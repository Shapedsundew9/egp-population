"""Gene pool management for Erasmus GP."""


from copy import deepcopy
from hashlib import blake2b
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
_SIG_NOT_NULL_FUNC = lambda x: x['signature'] is not None
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
_SIGNATURE_SQL = 'WHERE ({signature} in {matches}) IS TRUE'
_INITIAL_GC_SQL = ('WHERE {input_types} = {itypes} AND {inputs}={iidx} '
                      'AND {output_types} = {otypes} AND {outputs}={oidx} '
                      'AND ({exclude_column} NOT IN {exclusions}) IS NOT FALSE ORDER BY RANDOM() LIMIT {limit}')


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
    'expand_empty_tuple': True
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
        3. Persistence of updates to the gene pool.
        4. Local caching of GC's.

    The gene pool must be consistent i.e. no entry can depend on a genetic code
    that is not in the gene_pool.

    The primary difference with the genomic_library is the presence of transient data
    and the fast (and space saving) UID method of referencing active GC's. For the same
    reasons validation is not performed unless in debug mode.

    The gene_pool should be more local (faster to access) than the genomic library
    due to the heavy transaction load.

    The public member self.pool is the local cache of gGC's. It is a dictionary
    mapping both 'ref' and 'signature' kets to gGC's.
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
        self.pool = {}
        self._gl = genomic_library
        self._pool = table(config)
        self._update_str = UPDATE_STR.replace('__table__', config['table'])
        self._interface = None
        self.encode_value = self._pool.encode_value
        self.select = partial(self._gl.select, container='pkdict')
        if self._pool.raw.creator:
            self._pool.arbitrary_sql(sql_functions())


    def _reference(self, gcs):
        """Create local references for genomic library GCs.

        gcs is modified by this method.

        Args
        ----
        gcs (list(dict)): Genomic library recursive GC select results.

        Returns
        -------
        gcs
        """
        for gc in gcs:
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


    def _interface_hash(self, input_eps, output_eps):
        """Create a 64-bit hash of the population interface definition.

        Args
        ----
        input_eps (iterable(int)): Iterable of input EP types.
        output_eps (iterable(int)): Iterable of output EP types.

        Returns
        -------
        (int): 64 bit hash as a signed 64 bit int.
        """
        h = blake2b(digest_size=8)
        for i in input_eps:
            h.update(i.to_bytes(2, 'little'))
        for o in output_eps:
            h.update(o.to_bytes(2, 'little'))
        a = int.from_bytes(h.digest(), 'little')
        return (0x7FFFFFFFFFFFFFFF & a) - (a & (1 << 63))


    def initialize(self, inputs, outputs, exclusions=tuple(), num=1, vt=vtype.OBJECT):
        """Fetch or create num valid GC's with the specified inputs & outputs.

        The initial population is constructed
        GC's matching the criteria in the gene pool will be given preference over
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
        input_eps, xputs['itypes'], xputs['iidx'] = interface_definition(inputs, vt)
        output_eps, xputs['otypes'], xputs['oidx'] = interface_definition(outputs, vt)
        self._interface = self._interface_hash(input_eps, output_eps)

        # Find the GC's that match and then recursively pull them from the genomic library
        matches = tuple(m[0] for m in self._gl.library.raw.select(_INITIAL_GC_SQL, literals=xputs, columns=('signature',)))
        _logger.info(f'Found {len(matches)} GCs matching input and output criteria.')
        if _LOG_DEBUG:
            _logger.debug(f'GC signatures: {[sha256_to_str(m) for m in matches]}.')
        self.pull(matches)
        for signature in matches:
            self.pool[signature]['interface'] = self._interface
            self.pool[signature]['modified'] = True

        # If there was not enough fill the whole population create some new gGC's & mark them as individuals too.
        # This may require pulling new agc's from the genomic library through steady state exceptions
        # in the stabilise() in which case we need to pull in all dependents not already in the
        # gene pool.
        ggcs = []
        for _ in range(num - len(matches)):
            gc_list = stablise(self._gl, eGC(inputs=inputs, outputs=outputs, vt=vt))
            ggcs.append(gGC(gc_list[0], interface=self._interface, modified=True))
            ggcs.extend((gGC(gc, modified=True) for gc in gc_list[1:]))
        self.add(ggcs)

    def add(self, ggcs):
        """Add gGGs to the gene pool pulling any sub-GC's from the genomic library as needed.

        Args
        ----
        ggcs (iterable(gGC)): gGCs to be added. All referenced GCAs and GCBs must either exist
        in the ggcs iterable or the genomic library.
        """
        children = []
        for ggc in filter(_MODIFIED_FUNC, ggcs):
            self.pool[ggc['ref']] = ggc
            gca = ggc.get('gca', None)
            if gca is not None and gca not in self.pool:
                children.append(gca)
            gcb = ggc.get('gcb', None)
            if gcb is not None and gcb not in self.pool:
                children.append(gcb)
        self.pull(children)
        self.push()

    def pull(self, signatures):
        """Add aGCs and all sub-GC's recursively from the genomic library to the gene pool.

        aGC's are converted to gGC's.

        Args
        ----
        signatures (iterable(bytes[32])): Signatures to pull from the genomic library.
        """
        gcs = self._reference(self._gl.recursive_select(_SIGNATURE_SQL, literals={'matches': signatures}))
        self.pool = {gc['ref']: gGC(gc, modified=True) for gc in gcs}
        self.pool.update({gc['signature']: self.pool[gc['ref']] for gc in filter(_SIG_NOT_NULL_FUNC, self.pool.values())})
        self.push()

    def push(self):
        """Insert or update into locally modified gGC's into the persistent gene_pool."""
        # TODO: This can be optimised to further minimise the amount of data munging of unmodified values.
        modified_gcs = (gc for gc in filter(_MODIFIED_FUNC, self.pool.values()))
        for updated_gc in self._pool.upsert(modified_gcs, self._update_str, {}, _UPDATE_RETURNING_COLS):
            gc = self.pool[updated_gc['ref']]
            gc.update(updated_gc)
            for col in HIGHER_LAYER_COLS:
                gc[col] = gc[col[1:]]
        for gc in self.pool.values():
            gc['modified']= False

    def delete(self, refs):
        """Delete GCs from the gene pool.

        Args
        ----
        refs(iterable(int)): GC 'ref' values to delete.
        """
        for ref in refs:
            if self.pool[ref].get('signature', None) is not None:
                del self.pool[self.pool[ref]['signature']]
            del self.pool[ref]
        self._pool.delete('{ref} in {ref_tuple}', {'ref_tuple': tuple(refs)})