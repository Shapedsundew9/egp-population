from egp_types.reference import reference
from itertools import count
from copy import deepcopy
from functools import partial
from os.path import dirname, join
from pypgtable import table
from pypgtable.typing import TableSchema, Conversions, TableConfig
from json import load, loads, dumps
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger, Logger
from .typing import Populations, Population
from typing import Literal

_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
# TODO: Add a _LOG_CONSISTENCY which additionally does consistency checking
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)

_LAST_OWNER_ID_QUERY_SQL: str = 'WHERE uid = {uid}'
_LAST_OWNER_ID_UPDATE_SQL: str = "{last_owner_id} = (({last_owner_id}::BIGINT + 1::BIGINT) & x'7FFFFFFF::BIGINT)::INTEGER')"

_POPULATION_IN_DEFINITION: Literal[-1] = -1
_MAX_INITIAL_LIST: Literal[100000] = 100000
_MIN_PGC_LAYER_SIZE: Literal[100] = 100
_MAX_PGC_LAYER_SIZE: Literal[10000] = 10000
_MINIMUM_SUBPROCESS_TIME: Literal[60] = 60
_MINIMUM_AVAILABLE_MEMORY: Literal[134217728] = 128 * 1024 * 1024
_MAX_POPULATION_SIZE: Literal[100000] = 100000


with open(join(dirname(__file__), "formats/population_table_format.json"), "r") as file_ptr:
    _POPULATION_TABLE_SCHEMA: TableSchema = load(file_ptr)
with open(join(dirname(__file__), "formats/population_metrics_table_format.json"), "r") as file_ptr:
    _POPULATION_METRICS_TABLE_SCHEMA: TableSchema = load(file_ptr)
with open(join(dirname(__file__), "formats/spuid_table_format.json"), "r") as file_ptr:
    _SPUID_TABLE_SCHEMA: TableSchema = load(file_ptr)

_POPULATIONS_CONVERSIONS: Conversions = (
    ('inputs', dumps, loads),
    ('outputs', dumps, loads)
)


_DEFAULT_POPULATIONS_CONFIG: TableConfig = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool_populations',
    'schema': _POPULATION_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
    'conversions': _POPULATIONS_CONVERSIONS
}
_DEFAULT_GP_SPUID_CONFIG: TableConfig = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool_spuid',
    'schema': _SPUID_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
}
_DEFAULT_POPULATION_METRICS_CONFIG: TableConfig = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'population_metrics',
    'schema': _POPULATION_METRICS_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
}


class population():

    def __init(self) -> None:
        self._owner_counters: dict[int, count] = {}
        self.reference_function: partial[int] = partial(reference, counters=self._owner_counters)
        self.owner_id: int = self.next_owner_id()
        self._next_reference: partial[int] = partial(self.reference_function, owner=self.owner_id)
        self._populations_table: table = table(_DEFAULT_POPULATIONS_CONFIG)
        _logger.info('Population(s) established.')

    def next_owner_id(self) -> int:
        """Get the next available owner UID."""
        return next(self._populations_table.update(_LAST_OWNER_ID_UPDATE_SQL, _LAST_OWNER_ID_QUERY_SQL,
            returning=('last_owner_id',), container='tuple'))[0]

    def purge_owner_counters(self) -> None:
        """When the owners no longer exist clean up the counters"""
        self._owner_counters = {}

    def create_population(self, config:Dict = {}) -> Dict[str, Any]:
        """Create a population in the gene pool.

        The gene pool can support multiple populations.
        Once a population is created it cannot be modified.

        The population entry is initially created with a size of _POPULATION_IN_DEFINITION.
        This indicates the population is being created in the gene pool. Once complete the size
        is updated to the configured size.

        If the population entry already exists it will either be defined (gene pool populated)
        and the size entry will be set to something other than _POPULATION_IN_DEFINITION or
        it will be in definition in which case the method waits indefinately checking every
        1.0s (on average) until size changes from _POPULATION_IN_DEFINITION.

        Args
        ----
        config (dict): (required)
            'size':(int): Defines the population size.
            'name':(str): An unique short (<= 64 charaters) string giving a name.
            'inputs':iter(vtype): Target individual input definition using vt type.
            'outputs':iter(vtype): Target individual output definition using vt type.
            'characterize':callable(population): Characterize the population.
                The callable must take a single parameter that is an iterable of gGC's.
                The callable returns an iterable of characterisation dicts in the same order.
            'vt':(vtype): A vtype value defining how to interpret inputs & outputs.
            'description':(str): An arbitary string with a longer (<= 8192 characters) description.
            'meta_data':(str): Additional data about the population stored as a string (unlimited length).

        Returns
        -------
        data (dict):
            'idx':(int): A unique id given to the population.
            'worker_id:(int): The id of the worker that created the initial population.
            'size':(int): Defines the population size.
            'name':(str): An unique short (<= 64 charaters) string giving a name.
            'inputs':list(str): Target individual input definition as strings.
            'outputs':list(str): Target individual output definition as strings.
            'characterize':callable(gc): Characterize an individual of the population.
            'recharacterize': callable(g): Re-calculate survivability for a characterized individual.
            'vt':(vtype): vtype.EP_TYPE_STR.
            'description':(str): An arbitary string with a longer (<= 8192 characters) description.
            'meta_data':(str): Additional data about the population stored as a string (unlimited length).
            'created':(datetime): The time at which the population entry was created.
        """
        data = deepcopy(config)
        data['inputs'] = [asstr(i, config['vt']) for i in config['inputs']]
        data['outputs'] = [asstr(o, config['vt']) for o in config['outputs']]
        _, input_types, inputs = interface_definition(config['inputs'], config['vt'])
        _, output_types, outputs = interface_definition(config['outputs'], config['vt'])
        data['oih'] = ordered_interface_hash(input_types, output_types, inputs, outputs)
        data['size'] = _POPULATION_IN_DEFINITION
        data['worker_id'] = self.worker_id
        data = next(self._populations_table.insert((data,), '*'))
        data['vt'] = vtype.EP_TYPE_STR
        data['characterize'] = partial(_characterization_wrapper, characterize=config['characterize'])
        data['recharacterize'] = config['recharacterize']

        # If this worker created the entry then populate the gene pool.
        if data['worker_id'] == self.worker_id:
            _logger.info(f"Creating new population: {data['name']}.")
            self._population_data[data['uid']] = data
            self._populate(data['inputs'], data['outputs'], data['uid'], num=config['size'], vt=data['vt'])
            _logger.debug('Back to population create scope.')
            data['size'] = config['size']
            self._populations_table.update("{size} = {s}", "{uid} = {i}", {'s': data['size'], 'i': data['uid']})
        else:
            while data['size'] == _POPULATION_IN_DEFINITION:
                _logger.info(f"Waiting for population {data['name']} to be created by worker {data['worker_id']}.")
                sleep(0.5 + random())
                data = next(self._populations_table.select('{name} = {n}', {'n': data['name']}))
        return data

    def register_characterization(self, population_uid:int, characterize:Callable, recharacterize:Callable):
        """Register a function to characterize an individual & recharactize as the population changes.

        Every population in the Gene Pool that is to be evolved must have both a characterization
        function and a recharacterization function defined. The characterization function takes an individual
        gGC and calculates a fitness score and a survivability score. Both fitness
        & survivability are values between 0.0 and 1.0 inclusive. The recharacterization function
        calculates the survivability of a previously characterized gGG.

        Fitness is the fitness of the individual to the solution. A value of 1.0 means the
        solution is good enough & evolution ends.

        Survivability is often strongly corrolated with fitness but is not the same. Survivability
        is the relative weight of the indivdiual in the population for surviving to the next generation.
        It is different from fitness to allow diversity & novelty to be explored (which may lead
        unintuitively to greater fitness.)  The assumption is
        that survivability is a function of the population and not just the individual so as the
        population mutates the survivability of the individual changes.

        Args
        ----
        population_uid (int): The population UID to register the characterisation function to
        characterize (f(ggc)->tuple(float, float)
        recharacterize (f(ggc)->float)
        """
        self._population_data[population_uid]['characterize'] = characterize
        self._population_data[population_uid]['recharacterize'] = recharacterize

