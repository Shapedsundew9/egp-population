"""Gene pool management for Erasmus GP."""


from copy import deepcopy
from functools import partial
from os.path import dirname, join
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger
from json import load
from egp_genomic_library.genomic_library import compress_json, decompress_json, sql_functions, UPDATE_STR, UPDATE_RETURNING_COLS, HIGHER_LAYER_COLS, sha256_to_str
from egp_physics.ep_type import vtype
from egp_physics.gc_type import eGC, interface_definition, gGC, random_reference
from egp_physics.physics import stablise
from egp_physics.gc_graph import gc_graph

from pypgtable import table


_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)
_LOG_INFO = _logger.isEnabledFor(INFO)
_LOG_WARN = _logger.isEnabledFor(WARN)
_LOG_ERROR = _logger.isEnabledFor(ERROR)
_LOG_FATAL = _logger.isEnabledFor(FATAL)


_MODIFIED_FUNC = lambda x: x['modified']
_UPDATE_RETURNING_COLS = tuple((c for c in filter(lambda x: x != 'signature', UPDATE_RETURNING_COLS))) + ('ref',)


def compress_igraph(x):
    """Extract the internal representation and compress."""
    return compress_json(x.graph)


def decompress_igraph(x):
    """Create a gc_graph() from an decompressed internal representation."""
    return gc_graph(internal=decompress_json(x))


_CONVERSIONS = (
    ('graph', compress_json, decompress_json),
    ('meta_data', compress_json, decompress_json),
    ('igraph', compress_igraph, decompress_igraph)
)


_PTR_MAP = {
    'gca_ref': 'ref',
    'gcb_ref': 'ref'
}


with open(join(dirname(__file__), "formats/table_format.json"), "r") as file_ptr:
    _GP_TABLE_SCHEMA = load(file_ptr)


# GC queries
_SIGNATURE_SQL = 'WHERE {signature} in {matches}'
_INITIAL_GC_EX_SQL = ('WHERE {input_types}={itypes} AND {inputs}={iidx} '
                      'AND {output_types}={otypes} AND {outputs}={oidx} '
                      'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT {limit}')
_INITIAL_GC_SQL = ('WHERE {input_types}={itypes} AND {inputs}={iidx} '
                   'AND {output_types}={otypes} AND {outputs}={oidx} '
                   'ORDER BY RANDOM() LIMIT {limit}')


# The default config
_DEFAULT_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool',
    'ptr_map': _PTR_MAP,
    'schema': _GP_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
    'conversions': _CONVERSIONS,
}


def default_config():
    """Return a deepcopy of the default genomic library configuration.

    The copy may be modified and used to create a genomic library instance.

    Returns
    -------
    (dict): The default genomic_library() configuration.
    """
    return deepcopy(_DEFAULT_CONFIG)


class gene_pool():
    """Store of transient genetic codes & associated data for a population.

    The gene_pool is responsible for:
        1. Populating calculable entry fields.
        2. Providing an application interface to the fields.

    The gene pool must be consistent i.e. no entry can depend on a genetic code
    that is not in the gene_pool or genomic_library.

    The primary difference with the genomic_library is the presence of transient data
    and the fast (and space saving) UID method of referencing active GC's. For the same
    reasons validation is not performed unless in debug mode.

    The gene_pool should be more local (faster to access) than the genomic library
    due to the heavy transaction load.
    """

    def __init__(self, genomic_library, config=_DEFAULT_CONFIG):
        """Connect to or create a gene pool.

        The gene pool data persists in a postgresql database. Multiple
        instances of the gene_pool() class can connect to the same database
        (use the same configuration).

        Args
        ----
        genomic_library (genomic_library): Source of genetic material.
        config(pypgtable config): The config is deep copied by pypgtable.
        """
        self._gl = genomic_library
        self._pool = table(config)
        self._update_str = UPDATE_STR.replace('__table__', config['table'])
        self.encode_value = self._pool.encode_value
        self.select = partial(self._gl.select, container='pkdict')
        if self._pool.raw.creator:
            self._pool.arbitrary_sql(sql_functions())


    def _reference(self, gcs):
        """Create local references for genomic library GCs.

        gcs is modified by this method.

        Args
        ----
        gcs (pkdict): Genomic library recursive GC select results.

        Returns
        -------
        gcs
        """
        for gc in gcs.values():
            if 'ref' not in gc:
                gc['ref'] = random_reference()

            gca = gc['gca']
            if gca is not None and 'gca_ref' not in gc:
                if 'ref' not in gcs[gca]:
                    gcs[gca]['ref'] = random_reference()
                gc['gca_ref'] = gcs[gca]['ref']

            gcb = gc['gcb']
            if gcb is not None and 'gcb_ref' not in gc:
                if 'ref' not in gcs[gcb]:
                    gcs[gcb]['ref'] = random_reference()
                gc['gcb_ref'] = gcs[gcb]['ref']

        return gcs


    def initialize(self, inputs, outputs, exclusions=tuple(), num=1, vt=vtype.OBJECT):
        """Create num valid GC's with the specified inputs & outputs.

        GC's matching the criteria in the GMS will be given preference over
        creating new GC's. If there are more GC's in the GMS that meet the
        criteria than num then the returned GC's will be randomly selected.
        If there are less then valid GC's with the correct inputs and outputs
        will be created randomly using the GMS.

        If the input & output specification is such that no valid GC can be found
        or created None is returned.

        Args
        ----
        inputs (iterable(object or type)): Input objects of vt type.
        outputs (iteratble(object or type)): Output objects of vt type.
        gms (object): Genetic Material Source: A gene_pool or genomic_library object.
        num (int): The number of GC's to create
        vt (vtype): See vtype definition.

        Returns
        -------
        (list(gGC)) or None
        """
        xputs = {
            'exclude_column': 'signature',
            'exclusions': exclusions,
            'limit': num
        }
        _, xputs['itypes'], xputs['iidx'] = interface_definition(inputs, vt)
        _, xputs['otypes'], xputs['oidx'] = interface_definition(outputs, vt)

        # Find the GC's that match and then recursively pull them from the genomic library
        sql_str = _INITIAL_GC_EX_SQL if exclusions else _INITIAL_GC_SQL
        matches = tuple(m[0] for m in self._gl.library.raw.select(sql_str, literals=xputs, columns=('signature',)))
        _logger.info(f'Found {len(matches)} GCs matching input and output criteria.')
        if _LOG_DEBUG:
            _logger.debug(f'GC signatures: {[sha256_to_str(m) for m in matches]}.')
        if matches:
            gcs = self._reference(self._gl.select(_SIGNATURE_SQL, literals={'matches': matches}))
            ggcs = {gc['ref']: gGC(gc, individual=gc['signature'] in matches, modified=True) for gc in gcs.values()}
        else:
            ggcs = {}

        # If there was not enough fill the rest with gGC's & mark them as individuals too.
        for _ in range(num - len(matches)):
            ggc = gGC(stablise(self._gl, eGC(inputs=inputs, outputs=outputs, vt=vt)), individual=True, modified=True)
            ggcs[ggc['ref']] = ggc

        return self.upsert(ggcs)

    def upsert(self, gcs):
        """Insert or update into the gene_pool.

        This method modifies gcs.

        Args
        ----
        gcs (dict(dict)): keys are references and dicts are genetic code
            entries. Values will be normalised & updated in place.

        Returns
        -------
        gcs
        """
        # TODO: This can be optimised to further minimise the amount of data munging of unmodified values.
        modified_gcs = (gc for gc in filter(_MODIFIED_FUNC, gcs.values()))
        for updated_gc in self._pool.upsert(modified_gcs, self._update_str, {}, _UPDATE_RETURNING_COLS):
            gc = gcs[updated_gc['ref']]
            gc.update(updated_gc)
            gc['modified'] = False
            for col in HIGHER_LAYER_COLS:
                gc[col] = gc[col[1:]]
        return gcs

    def delete(self, refs):
        """Delete GCs from the gene pool.

        Args
        ----
        refs(iterable(int)): GC 'ref' values to delete.
        """
        self._pool.delete('{ref} in {ref_tuple}', {'ref_tuple': tuple(refs)})