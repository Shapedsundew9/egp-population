"""Gene pool management for Erasmus GP."""


from copy import deepcopy
from os.path import dirname, join
from logging import DEBUG, NullHandler, getLogger
from json import load
from egp_genomic_library.genomic_library import compress_json, decompress_json, sql_functions
from egp_physics.ep_type import vtype
from egp_physics.gc_type import eGC, interface_definition, gGC
from pypgtable import table


_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOGIT = _logger.getEffectiveLevel() == DEBUG


_CONVERSIONS = (
    ('graph', compress_json, decompress_json),
    ('meta_data', compress_json, decompress_json),
    ('igraph', compress_json, decompress_json)
)


_PTR_MAP = {
    'gca_ref': 'ref',
    'gcb_ref': 'ref'
}


with open(join(dirname(__file__), "formats/table_format.json"), "r") as file_ptr:
    _GP_TABLE_SCHEMA = load(file_ptr)
_HIGHER_LAYER_COLS = tuple((key for key in filter(lambda x: x[0] == '_', _GP_TABLE_SCHEMA)))
_UPDATE_RETURNING_COLS = tuple((x[1:] for x in _HIGHER_LAYER_COLS)) + ('updated', 'created')


# Initial GC query
_INITIAL_GC_SQL = ('WHERE "input_types"={input_types} AND "inputs"={inputs} '
                   'AND "output_types"={output_types} AND "outputs"={outputs} '
                   'AND {exclude_column} NOT IN {exclusions} ORDER BY RANDOM() LIMIT {limit}')


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
        if self._pool.raw.creator:
            with open(join(dirname(__file__), 'data/gl_functions.sql'), 'r') as fileptr:
                self._pool.arbitrary_sql(sql_functions())

    def encode_value(self, column, value):
        """Encode value using the registered conversion function for column.

        If no conversion function is registered value is returned unchanged.
        This function is provided to create encoded literals in query functions.

        Args
        ----
        column (str): Column name for value.
        value  (obj): Value to encode.

        Returns
        -------
        (obj): Encoded value
        """
        return self._pool.encode_value(column, value)

    def initialize(self, inputs, outputs, exclusions=tuple(), num=1, vt=vtype.OBJECT):
        """Create num valid GC's with the specified inputs & outputs.

        GC's matching the criteria in the GMS will be given preference over
        creating new GC's. If there are more GC's in the GMS that meet the
        criteria than num then the returned GC's will be randomly selected.
        If there are less then valid GC's with the correct inputs and outputs
        will be created randomly.

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
        (list(eGC)) or None
        """
        xputs = {
            'exclude_column': 'signature',
            'exclusions': exclusions,
            'limit': num
        }
        _, xputs['input_types'], xputs['inputs'] = interface_definition(inputs, vt)
        _, xputs['output_types'], xputs['outputs'] = interface_definition(outputs, vt)

        # Find the GC's that match and then recursively pull them from the genomic library
        matches = self._gl.library.raw.select(_INITIAL_GC_SQL, literals=xputs, columns=('signature',))
        _logger.info(f'Found {len(matches)} GCs matching input and output criteria.')
        if _LOGIT:
            _logger.debug(f'GC signatures: {matches}.')
        gcs = self._gl.select('WHERE signature in {matches}', literals={'matches': matches})

        # All GCs pulled from the genomic library need thier higher layer fields set.
        for gc in gcs.values():
            for hlk in _HIGHER_LAYER_COLS:
                gc[hlk] = gc[hlk[1:]]

        # GC's that matched are individuals - the rest genetic flotsam and jetsam
        for signature in matches:
            gcs[signature]['individual'] = True

        # If there was not enough fill the rest with gGC's & mark them as individuals too.
        new_gcs = (eGC(inputs=inputs, outputs=outputs, vt=vt) for _ in range(num - len(matches)))
        # TODO: ensure eGC is in a steady state then make it a gGC with the individual flag set.

        gcs.update({})
        return gcs

    def select(self, query_str='', literals={}):
        """Select genetic codes from the gene pool.

        The select is done recursively returning all constituent genetic codes and
        codons.

        Args
        ----
        query_str (str): Query SQL: e.g. '{input_types} @> {inputs}'.
        literals (dict): Keys are labels used in query_str. Values are literals to replace the labels.
            NOTE: Literal values for encoded columns must be encoded. See encode_value().

        Returns
        -------
        (dict(dict)): 'signature' keys with gc values.
        """
        return self._library.recursive_select(query_str, literals, container='pkdict')

    def upsert(self, entries):
        """Insert or update into the gene_pool.

        Args
        ----
        entries (dict(dict)): keys are references and dicts are genetic code
            entries. Values will be normalised & updated in place
        """
        self._normalize(entries)
        for_insert = (x for x in filter(lambda x: not x.get('_stored', False), entries.values()))
        for_update = (x for x in filter(lambda x: x.get('_stored', False), entries.values()))
        updated_entries = self._library.upsert(for_insert, _UPDATE_STR, {}, _UPDATE_RETURNING_COLS, _HIGHER_LAYER_COLS)
        updated_entries.extend(self._library.update(for_update, _UPDATE_STR, {}, _UPDATE_RETURNING_COLS))
        for updated_entry in updated_entries:
            entry = entries[updated_entry['signature']]
            entry.update(updated_entry)
            for hlk in _HIGHER_LAYER_COLS:
                entry[hlk] = entry[hlk[1:]]
        for entry in entries:
            entry['_stored'] = True
