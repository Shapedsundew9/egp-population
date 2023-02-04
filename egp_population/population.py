from egp_types.reference import reference
from itertools import count
from copy import deepcopy
from functools import partial
from os.path import dirname, join
from pypgtable import table
from pypgtable.typing import TableSchema, Conversions, TableConfigNorm, TableConfig
from pypgtable.validators import raw_table_config_validator
from json import load, loads, dumps
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger, Logger
from .typing import PopulationsConfig, PopulationsConfigNorm, PopulationConfig, PopulationConfigNorm
from typing import Literal
from .population_validator import population_entry_validator


_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
# TODO: Add a _LOG_CONSISTENCY which additionally does consistency checking
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


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


_POPULATIONS_CONVERSIONS: Conversions = (
    ('inputs', dumps, loads),
    ('outputs', dumps, loads)
)


_DEFAULT_POPULATIONS_CONFIG: TableConfigNorm = raw_table_config_validator.normalized({
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool_populations',
    'schema': _POPULATION_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
    'conversions': _POPULATIONS_CONVERSIONS
})

_DEFAULT_POPULATION_METRICS_CONFIG: TableConfigNorm = raw_table_config_validator.normalized({
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'population_metrics',
    'schema': _POPULATION_METRICS_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
})


class population():

    def __init__(self, populations_config: PopulationsConfig | PopulationsConfigNorm,
        table_config: TableConfig | TableConfigNorm = _DEFAULT_POPULATIONS_CONFIG) -> None:
        
        self.table: table = table(raw_table_config_validator.normalized(table_config))
        self.configs: dict[int, PopulationConfigNorm] = {p['uid']: population_entry_validator.normalized(p) for p in populations_config}
        _logger.info('Population(s) established.')

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

    def _new_population(self, population: PopulationNorm, num: int) -> None:
        """Fetch or create num target GC's & recursively pull in the pGC's that made them.

        Construct a population with inputs and outputs as defined by inputs, outputs & vt.
        Inputs & outputs define the GC's interface.

        There are 2 sources of valid individuals: The higher layer (GL in this case) or to
        create them.
        FIXME: Add pull from higher layer (recursively).
        FIXME: Implement forever work.

        Args
        ----
        population: Population definition.
        num: The number of GC's to create
        vt: See vtype definition.
        """
        # Check the population index exists
        _logger.info(f"Adding {num} GC's to population index: {population['uid']}.")

        # If there was not enough fill the whole population create some new gGC's & mark them as individuals too.
        # This may require pulling new agc's from the genomic library through steady state exceptions
        # in the stabilise() in which case we need to pull in all dependents not already in the
        # gene pool.
        _logger.info(f'{num} GGCs to create.')
        for _ in range(num):
            egc: eGC = eGC(inputs=population['inputs'], outputs=population['outputs'], vt=vtype.EP_TYPE_STR)
            rgc, fgc_dict = stablize(self._gl, egc)

            # Just in case it is trickier than expected.
            retry_count: int = 0
            while rgc is None:
                _logger.info('eGC random creation failed. Retrying...')
                retry_count += 1
                egc = eGC(inputs=population['inputs'], outputs=population['outputs'], vt=vtype.EP_TYPE_STR)
                rgc, fgc_dict = stablize(self._gl, egc)
                if retry_count == 3:
                    raise ValueError(f"Failed to create eGC with inputs = {population['inputs']} and outputs"
                                     f" = {population['outputs']} {retry_count} times in a row.")

            rgc['population_uid'] = population['uid']
            self.pool[rgc['ref']] = rgc
            self.pool.update(fgc_dict)
            _logger.debug(f'Created GGCs to add to Gene Pool: {[ref_str(ggc["ref"]) for ggc in (rgc, *fgc_dict.values())]}')cd ../egp   