"""Gene pool management for Erasmus GP."""

from copy import deepcopy
from functools import partial
from os.path import dirname, join
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger
from json import load, loads, dumps
from random import random
from numpy.core.fromnumeric import mean
from numpy.random import choice
from numpy import array, float32, sum
from egp_genomic_library.genetic_material_store import genetic_material_store
from egp_genomic_library.genomic_library import (compress_json, decompress_json, UPDATE_STR, UPDATE_RETURNING_COLS,
    sql_functions, HIGHER_LAYER_COLS)
from egp_physics.ep_type import asstr, vtype, asint
from egp_physics.gc_type import eGC, interface_definition, gGC, interface_hash
from egp_physics.physics import cull_physical, stablise, xGC_evolvability, pGC_fitness, select_pGC, xGC_inherit, random_reference, RANDOM_PGC_SIGNATURE
from egp_physics.gc_graph import gc_graph
from egp_physics.execution import remove_callable
from .gene_pool_cache import gene_pool_cache, SP_UID_LOWER_LIMIT, SP_UID_UPPER_LIMIT
from time import time, sleep
from multiprocessing import set_start_method, Process
from os import cpu_count, kill
from psutil import Process as psProcess
from psutil import virtual_memory
from gc import enable, disable, collect, freeze, unfreeze
from signal import signal, SIGUSR1

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
_POPULATION_IN_DEFINITION = -1
_MAX_INITIAL_LIST = 10000
_MIN_PGC_LAYER_SIZE = 100
_MINIMUM_SUBPROCESS_TIME = 60
_MINIMUM_AVAILABLE_MEMORY = 128 * 1024 * 1024


def compress_igraph(x):
    """Extract the internal representation and compress."""
    return compress_json(x.graph)


def decompress_igraph(x):
    """Create a gc_graph() from an decompressed internal representation."""
    return gc_graph(internal=decompress_json(x))


_GP_CONVERSIONS = (
    ('graph', compress_json, decompress_json),
    ('meta_data', compress_json, decompress_json),
    ('igraph', compress_igraph, decompress_igraph)
)
_POPULATIONS_CONVERSIONS = (
    ('inputs', dumps, loads),
    ('outputs', dumps, loads)
)


# Tree structure
_LEL = 'gca_ref'
_REL = 'gcb_ref'
_NL = 'ref'
_PTR_MAP = {
    _LEL: _NL,
    _REL: _NL
}


with open(join(dirname(__file__), "formats/gp_table_format.json"), "r") as file_ptr:
    _GP_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/gp_spuid_table_format.json"), "r") as file_ptr:
    _GP_SPUID_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/gpp_table_format.json"), "r") as file_ptr:
    _GPP_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/gp_metrics_format.json"), "r") as file_ptr:
    _GP_METRICS_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/layer_metrics_format.json"), "r") as file_ptr:
    _LAYER_METRICS_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/pgc_metrics_format.json"), "r") as file_ptr:
    _PGC_METRICS_TABLE_SCHEMA = load(file_ptr)


# GC queries
_INITIAL_GP_COLUMNS = ('ref', 'evolvability', 'signature')
_INITIAL_GL_COLUMNS = ('signature', 'evolvability')
_SP_UID_UPDATE_STR = '{next_spuid} = COALESCE(NULLIF({next_spuid} + 1, {max_spuid}), {next_spuid})'
_SIGNATURE_SQL = 'WHERE ({signature} = ANY({matches}))'
_INITIAL_GC_SQL = ('WHERE {input_types} = {itypes} AND {inputs}={iidx} '
                      'AND {output_types} = {otypes} AND {outputs}={oidx} '
                      'AND NOT({exclude_column} = ANY({exclusions})) ORDER BY RANDOM() LIMIT {limit}')
_PGC_DEFINITION = {
    'itypes': [asint('egp_physics.gc_type_gGC')],
    'otypes': [asint('egp_physics.gc_type_gGC')],
    'iidx': b'0',
    'oidx': b'0',
    'exclude_column': 'signature',
    'exclusions': [RANDOM_PGC_SIGNATURE],
    'limit': _MAX_INITIAL_LIST
}
_COMPETITION_SQL = 'WHERE {population} = {popidx} and {ref} not in {exclusions}'


# The gene pool default config
_DEFAULT_GP_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool',
    'ptr_map': _PTR_MAP,
    'schema': _GP_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
    'conversions': _GP_CONVERSIONS
}
_DEFAULT_POPULATIONS_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool_populations',
    'schema': _GPP_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
    'conversions': _POPULATIONS_CONVERSIONS
}
_DEFAULT_GP_SPUID_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gene_pool_spuid',
    'schema': _GP_SPUID_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
}
_DEFAULT_GP_METRICS_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'gp_metrics',
    'schema': _GP_METRICS_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
}
_DEFAULT_LAYER_METRICS_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'layer_metrics',
    'schema': _LAYER_METRICS_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
}
_DEFAULT_PGC_METRICS_CONFIG = {
    'database': {
        'dbname': 'erasmus'
    },
    'table': 'pgc_metrics',
    'schema': _PGC_METRICS_TABLE_SCHEMA,
    'create_table': True,
    'create_db': True,
}
_DEFAULT_CONFIGS = {
    "gp": _DEFAULT_GP_CONFIG,
    "populations": _DEFAULT_POPULATIONS_CONFIG,
    "gp_spuid": _DEFAULT_GP_SPUID_CONFIG,
    "gp_metrics": _DEFAULT_GP_METRICS_CONFIG,
    "layer_metrics": _DEFAULT_LAYER_METRICS_CONFIG,
    "pgc_metrics": _DEFAULT_PGC_METRICS_CONFIG,
}


def default_config():
    """Return a deepcopy of the default genomic library configuration.

    The copy may be modified and used to create a genomic library instance.

    Returns
    -------
    (dict): The default genomic_library() configuration.
    """
    return deepcopy(_DEFAULT_CONFIGS)


class gene_pool(genetic_material_store):
    """Store of transient genetic codes & associated data for a population.

    The gene_pool is responsible for:
        1. Populating calcuble entry fields.
        2. Providing an application interface to the fields.
        3. Persistence of updates to the gene pool.
        4. Local caching of GC's.

    The gene pool must be consistent i.e. no entry can depend on a genetic code
    that is not in the gene_pool.

    The primary difference with the genomic_library is the presence of transient data
    and the fast (and space saving) UID method of referencing and storing active GC's.
    For the same reasons validation is not performed unless in debug mode.

    The gene_pool is more local (faster to access) than the genomic library
    due to the heavy transaction load. It is also designed to be multi-process
    memory efficient.

    The public member self.pool is the local cache of gGC's.
    """

    #TODO: Default genomic_library should be local host
    def __init__(self, genomic_library, configs=_DEFAULT_CONFIGS):
        """Connect to or create a gene pool.

        The gene pool data persists in a postgresql database. Multiple
        instances of the gene_pool() class can connect to the same database
        (use the same configuration).

        Args
        ----
        genomic_library (genomic_library): Source of genetic material.
        config(pypgtable config): The config is deep copied by pypgtable.
        """
        super().__init__(node_label=_NL, left_edge_label=_LEL, right_edge_label=_REL)

        # pool is a dictionary-like object of ref:gGC
        # All gGC's in pool must be valid & any modified are sync'd to the gene pool table
        # in the database at the end of the target epoch.
        self.pool = gene_pool_cache()
        _logger.info('Created Gene Pool Cache.')
        self._gl = genomic_library
        self._pool = table(configs['gp'])
        _logger.info('Established connection to Gene Pool table.')

        # UID for this worker
        self.worker_id = random_reference()
        _logger.info(f'Worker ID: {self.worker_id:08X}')

        # FIXME: Validators for all tables.

        # A dictionary of metric tables. Metric tables are undated at the end of the
        # taregt epoch.
        self._metrics = {c: table(configs[c]) for c in configs.keys() if 'metrics' in c}
        _logger.info('Established connections to metric tables.')

        # ?
        self._env_config = None
        self._population_data = {}
        self._populations_table = table(configs['populations'])
        if self._populations_table.raw.creator:
            _logger.info(f'This worker ({self.worker_id:08X}) is the Populations table creator.')
        self.populations()
        _logger.info('Populations(s) established.')

        # Modify the update strings to use the right table for the gene pool.
        self._update_str = UPDATE_STR.replace('__table__', configs['gp']['table'])

        #?
        self._interface = None
        self.encode_value = self._pool.encode_value
        self.select = self._pool.select

        # The number of layers in the gene pool
        # Layer 0 is the target layer
        # Layers 1 to N are pGC layers
        self.num_layers = 0
        self.pgc_layer_size = _MIN_PGC_LAYER_SIZE

        # Used to track the number of updates to individuals in a layer (?)
        self.layer_evolutions = [0]

        # If this instance created the gene pool then it is responsible for configuring
        # setting up database functions and initial population.
        if self._pool.raw.creator:
            _logger.info(f'This worker ({self.worker_id:08X}) is the Gene Pool table creator.')
            self._pool.raw.arbitrary_sql(sql_functions(), read=False)

        # The Sub-Process UID table
        self._spuid_table = table(configs['gp_spuid'])
        if self._spuid_table.raw.creator:
            self._spuid_table.insert(({'next_spuid': SP_UID_LOWER_LIMIT},))

        # The self terminate flag.
        # Set by the SIGUSR1 handler to True allowing the sub-process to exit gracefully writing its data
        # to the Gene Pool table when asked.
        self._terminate = False
        signal(SIGUSR1, self.self_terminate)


    def self_terminate(self):
        """Set the self termination flag.

        This is called by the SIGTERM handler.
        """
        self._terminate = True


    def get_spuid(self):
        """Get the next Sub-Process UID.

        The SPUID value must be atomically incremented for the next sub-process.
        """
        result = self._spuid_table.update(_SP_UID_UPDATE_STR, "*", {'max_spuid': SP_UID_UPPER_LIMIT}, returning="*")
        spuid = next(result)["next_spuid"] - 1
        _logger.info(f'Sub-process UID claimed: {spuid:08X}')
        return spuid


    def evolve(self, configs, num_sub_processes=None):
        """Co-evolve the population in pop_list.

        An evolution configuration = {
            'idx': Index of the population to evolve.
            'max_duration': Maximum duration (wall clock time) to evolve the population in seconds.
            'max_fitness': Stop evolution when at least one individual meets or exceeds this fitness.
            'max_epochs': Stop evolution when at least this number of epochs have passed.
        }

        Evolution of a population will stop when at least one of the criteria are met.
        This function will return only when all populations have met at least one of the critieria.

        Args
        ----
        configs (list(dict)): List of evolution configurations for co-evolved populations.
        num_sub_processes (int): Number of sub processes to spawn. If None the number of CPUs-1
                                 will be used.

        Returns
        -------
        (list(gGC)): List of fittest individual in each population.
        """
        self.pre_evolution_checks()
        while not self._exit_criteria():
            self._spawn(configs, num_sub_processes)
            self._purge_local_cache()
            self._repopulate_local_cache()


    def _spawn(self, num_sub_processes=None):
        """Spawn subprocesses.

        Args
        ----
        num_sub_processes (int): Number of sub processes to spawn. If None the number of CPUs-1
                                 will be used.
        """
        self._disconnect_tables()
        set_start_method("fork")
        if num_sub_processes is None:
            num_sub_processes = cpu_count() - 1

        # Memory clean up & freeze
        collect()
        disable()
        freeze()

        processes = [Process(target=self._entry_point) for _ in range(num_sub_processes)]
        start = time()
        for p in processes:

            # Sub-processes exit if the parent exits.
            p.daemon = True
            p.start()

        # Wait for a sub-process to terminate.
        while all((p.is_alive() for p in processes)) and self._memory_ok(start):
            sleep(1)

        # At least one sub-process is done or we need to reclaim some memory
        # so ask any running processes to finish up
        for p in filter(lambda x: x.is_alive(), processes):
            kill(p.pid, SIGUSR1)
        for p in processes:
            p.join()

        # Re-enable GC
        unfreeze()
        enable()

        self._reconnect_tables()

        # TODO: Are we done? Did we run out of sub-process IDs?


    def _memory_ok(self, start):
        """"Check if, after a minimum runtime, memory available is low.

        Only if we have been
        running for more than _MINIMUM_SUBPROCESS_TIME seconds do we check to see
        if RAM available is low (< _MINIMUM_AVAILABLE_MEMORY bytes). If it is
        low return False else True.

        Args
        ----
        start (float): The epoch time the sub-processes were spawned.

        Returns
        -------
        (bool).
        """
        duration = int(time() - start)
        if duration < _MINIMUM_SUBPROCESS_TIME:
            return True
        available = virtual_memory().available
        ok = available > _MINIMUM_AVAILABLE_MEMORY
        if not ok:
            _logger.info(f'Available memory is low ({available} bytes) after {duration}s.'
                          ' Signalling sub-processes to terminate.')
        return ok


    def _entry_point(self):
        """Entry point for sub-processes."""
        _logger.info("New sub-process created.")
        self._reconnect_tables()
        self.pool.set_sp_uid(self.get_spuid())
        if self.pool.sp_uid == SP_UID_UPPER_LIMIT - 2:
            _logger.info('Ran out of sub-process UIDs. Gene Pool must be purged and recreated to continue.')
            self.self_terminate()
        while not self._terminate:
            self.epoch()
        self._disconnect_tables()
        _logger.info("Sub-process gracefully terminating.")


    def populations(self):
        """Return the definition of all populations.

        The gene pool object stores a local copy of the population data with
        fitness & diversity functions defined (in 0 to all) cases.
        """
        self._population_data.update({p['idx']: p for p in self._populations_table.select()})
        if _LOG_DEBUG:
            _logger.debug('Populations table:\n'+str("\n".join(self._population_data.values())))
        return {p['name']: p for p in self._population_data}

    def create_population(self, config={}):
        """Create a population in the gene pool.

        The gene pool can support multiple populations.
        Once a population is created it cannot be modified.

        The population entry is initially created with a size of _POPULATION_IN_DEFINITION.
        This indicates the population is being created in the gene pool. Once complete the size
        is updated to the configured size.

        If the population entry already exists it will either be defined (gene pool populated)
        and the size entry will be set to something other than _POPULATION_IN_DEFINITION or
        it will be in definition in which case the method waits indefinately checking every
        1.5s (on average) until size changes from _POPULATION_IN_DEFINITION.

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
            'characterize':callable(population): Characterize the population.
            'vt':(vtype): vtype.EP_TYPE_STR.
            'description':(str): An arbitary string with a longer (<= 8192 characters) description.
            'meta_data':(str): Additional data about the population stored as a string (unlimited length).
            'created':(datetime): The time at which the population entry was created.
        """
        data = deepcopy(config)
        data['inputs'] = [asstr(i, config['vt']) for i in config['inputs']]
        data['outputs'] = [asstr(o, config['vt']) for o in config['outputs']]
        data['size'] = _POPULATION_IN_DEFINITION
        data['worker_id'] = self.worker_id
        data = next(self._populations_table.insert((data,), '*'))
        data['vt'] = vtype.EP_TYPE_STR

        # If this worker created the entry then populate the gene pool.
        if data['worker_id'] == self.worker_id:
            _logger.info(f"Creating new population: {data['name']}.")
            self._populate(data['inputs'], data['outputs'], data['idx'], num=config['size'], vt=data['vt'])
            data['size'] = config['size']
            self._populations_table.update("{size} = {s}", "{idx} = {i}", {'s': data['size'], 'i': data['idx']})
        else:
            while data['size'] == _POPULATION_IN_DEFINITION:
                _logger.info(f"Waiting for population {data['name']} to be created by worker {data['worker_id']}.")
                sleep(0.5 + random())
                data = next(self._populations_table.select('{name} = {n}', {'n': data['name']}))
        return data

    def register_characterization(self, pop_idx, characterization_func):
        """Register a function to characterize a population.

        Every population in the Gene Pool that is to be evolved must have a characterization
        function defined. The characterization function takes a population (iterable of
        individual gGC's) and calculates a fitness score and a survivability score. Both fitness
        & survivability are values between 0.0 and 1.0 inclusive.

        Fitness is the fitness of the individual to the solution. A value of 1.0 means the
        solution is good enough & evolution ends.

        Survivability is often strongly corrolated with fitness but is not the same. Survivability
        is the relative weight of the indivdiual in the population for surviving to the next epoch.
        It is different from fitness to allow diversity & novelty to be explored (which may lead
        unintuitively to greater fitness.)





        """

    def _populate(self, inputs, outputs, pop_idx, exclusions=[], num=1, vt=vtype.OBJECT):
        """Fetch or create num target GC's & recursively pull in the pGC's that made them.

        Construct a population with inputs and outputs as defined by inputs, outputs & vt.
        Inputs & outputs define the GC's interface.

        Selection of the target population:
            1. Find all GC's with the target interface in the GP & GL (to _MAX_INITIAL_LIST each)
            2. Find the maximum evolvability of each GC i.e. max at any layer.
            3. If > 2*num GC's in total: Randomly select num GC's weighted by max evolvability (no duplicates).
            4. If < num GC's in total: Add num - total randomly generated GC's.
            5. Characterise the target GC's.
            6. Cull to population size.

        Selection of pGC's
            1. Set the pGC layer size to the maximum population size (of all populations) or
               _MIN_PGC_LAYER_SIZE whichever is greater.
            2. If there are no pGCs in the GP
                a. Add the random pGC from the GL (this defines the initial self.num_layers)
                b. Randomly select (weighted by fitness in layer) layer size - 1 pGC's.
                   If there are <= layer size - 1 candidates all are added.

        Args
        ----
        inputs (iterable(object or type)): Input objects of vt type.
        outputs (iteratble(object or type)): Output objects of vt type.
        pop_idx (int): The population index to use for the new individuals.
        exclusions (iter(sha256)): GC signatures to exclude.
        num (int): The number of GC's to create
        vt (vtype): See vtype definition.
        """
        # Check the population index exists
        _logger.info(f"Adding {num} GC's to population index: {pop_idx}.")

        # Create search criteria for the genomic library
        xputs = {
            'exclude_column': 'signature',
            'exclusions': exclusions,
            'limit': _MAX_INITIAL_LIST
        }
        input_eps, xputs['itypes'], xputs['iidx'] = interface_definition(inputs, vt)
        output_eps, xputs['otypes'], xputs['oidx'] = interface_definition(outputs, vt)
        self._interface = interface_hash(input_eps, output_eps)

        # Find the GC's that match in the GP & GL
        matches = list(self._pool.raw.select(_INITIAL_GC_SQL, xputs, _INITIAL_GP_COLUMNS))
        gp_num = len(matches)
        _logger.info(f'Found {gp_num} GCs matching input and output criteria in the Gene Pool table.')
        matches.extend(self._gl.library.raw.select(_INITIAL_GC_SQL, xputs, _INITIAL_GL_COLUMNS))
        _logger.info(f'Found {len(matches) - gp_num} GCs matching input and output criteria in the genomic library.')

        # Matches is a list of (ref or signature, evolvability, )
        selected = [match[0] for match in matches]
        if len(matches) > 2*num:
            weights = array([max(match[1] for match in matches)], float32)
            weights /= sum(weights)
            selected = choice(selected, 2*num, False, weights)

        # For those selected from the GP mark them as individuals of this population
        # If they are already part of a population make a copy and give it a new reference.
        for uid in filter(lambda x: isinstance(x, int), selected):
            if self.pool[uid]['population']:
                clone = deepcopy(self.pool[uid])
                uid = self.pool.ref()
                self.pool[uid] = clone
            clone['population'] = pop_idx

        # For those from the genomic library pull them in and mark them as modified (new)
        # to the gene pool.
        from_gl = [uid for uid in filter(lambda x: not isinstance(x, int), selected)]
        self.pull_from_gl(from_gl)
        for signature in from_gl:
            self.pool[signature]['modified'] = True
            self.pool[signature]['population'] = pop_idx

        # GC's from the GL

        # If there was not enough fill the whole population create some new gGC's & mark them as individuals too.
        # This may require pulling new agc's from the genomic library through steady state exceptions
        # in the stabilise() in which case we need to pull in all dependents not already in the
        # gene pool.
        ggcs = []
        _logger.info(f'{num - len(matches)} GGCs to create.')
        for _ in range(num - len(matches)):
            rgc, fgc_dict = stablise(self._gl, eGC(inputs=inputs, outputs=outputs, vt=vt))
            ggcs.append(gGC(rgc, interface=self._interface, modified=True, population=pop_idx))
            ggcs.extend((gGC(gc, modified=True) for gc in fgc_dict.values()))
        _logger.debug(f'Created GGCs to add to Gene Pool: {[ggc["ref"] for ggc in ggcs]}')
        self.add_to_gp_cache(ggcs)
        self.push_to_gp()

        # Make an initial assessment of fitness
        for individual in self.individuals(pop_idx):
            individual['fitness'][0] = self._population_data[pop_idx]['fitness'](individual)
            individual['f_count'][0] = 1

        # Purge any overstock
        self.cull_target(pop_idx)

        # Populate pGC's
        self.pgc_layer_size = max((self.pgc_layer_size, max((p['size'] for p in self._population_data.values()))))
        pgcs = [pgc['ref'] for pgc in self._pool.values() if pgc['population'] == 0]
        if not len(pgcs):
            self.pull_from_gl((RANDOM_PGC_SIGNATURE,))
            self.num_layers = len(self.pool[RANDOM_PGC_SIGNATURE]['fitness'])
            pgcs = list(self._gl.select(_INITIAL_GC_SQL, _PGC_DEFINITION, ('signature', 'f_count', 'fitness')))

            # Extend the layers for all pGC's up to the number of layer in the random pGC.
            for pgc in pgcs:
                if len(pgc['f_count']) < self.num_layers:
                    pgc['f_count'].extend([0] * (self.num_layers - len(pgc['f_count'])))
                    pgc['fitness'].extend([0.0] * (self.num_layers - len(pgc['f_count'])))
                if len(pgc['f_count']) > self.num_layers:
                    _logger.fatal(f"Consistency error: pGC {pgc['signature']} has more layers than the random pGC.")
                    assert len(pgc['f_count']) > self.num_layers

            # Use the standard cull function to identify the excess.
            # Only pGC's that are used in at least one layer are pulled from the GL.
            for layer in range(1, self.num_layers):
                cull_physical(pgcs, layer)
            self.pull_from_gl()(pgc['signature'] for pgc in pgcs if any(pgc['f_valid']))

    def add_to_gp_cache(self, ggcs):
        """Add gGGs to the local gene pool cache pulling any sub-GC's from the genomic library as needed.

        NOTE: This *MUST* be the only function pushing GC's to the local GP cache (aside from self.pull_from_gl())
        NOTE: This function does *NOT* update the persistent GP. The intention is to do that (slow) operation infrequently.

        Args
        ----
        ggcs (iterable(gGC)): gGCs to be added. All referenced GCAs and GCBs must either exist
        in the ggcs iterable, gp or the genomic library (i.e. have a signature).
        """
        children = []
        for ggc in ggcs:
            # Muchos consistency checking.
            # This function is critical as we need to ensure:
            #   a) Everything going into the local GP cache meets a minimum bar.
            #   b) It is efficient
            if _LOG_DEBUG:
                _logger.debug(f'Adding {ggc["ref"]} to local cached gene pool.')
                if not isinstance(ggc, gGC):
                    _logger.error(f"GC ref {ggc['ref']} is not a gGC it is a {type(ggc)}!")
                if ggc['ref'] in self.pool:
                    _logger.warning(f"GC ref {ggc['ref']} already in local GP cache. No-op but inefficient.")
                gca_ref = ggc.get('gca_ref', None)
                if gca_ref not in ggcs and gca_ref not in self.pool and ggc.get('gca', None) is None:
                    _logger.fatal(f"Consistency failure. gca_ref in ggc: {ggc} does not exist!")
                gcb_ref = ggc.get('gcb_ref', None)
                if gcb_ref not in ggcs and gcb_ref not in self.pool and ggc.get('gca', None) is None:
                    _logger.fatal(f"Consistency failure. gcb_ref in ggc: {ggc} does not exist!")

            # Assumptions:
            #   a) Is a gGC
            #   b) Is not in the GP
            #   c) If signature, gca or gcb is defined equivilent _ref is defined.
            ggc['modified'] = True
            self.pool[ggc['ref']] = ggc
            gca = ggc.get('gca', None)
            if gca is not None and ggc['gca_ref'] not in self.pool:
                children.append(gca)
                if _LOG_DEBUG: _logger.debug(f'Appended {gca} to list to pull from the Genomic Library.')
            gcb = ggc.get('gcb', None)
            if gcb is not None and ggc['gcb_ref'] not in self.pool:
                children.append(gcb)
                if _LOG_DEBUG: _logger.debug(f'Appended {gcb} to list to pull from the Genomic Library.')

        # Must pull dependents from GL before adding nodes to the graph.
        # add_nodes() assumes the graph is complete.
        self.pull_from_gl(children)
        self.add_nodes(ggcs)

    def pull_from_gl(self, signatures):
        """Pull aGCs and all sub-GC's recursively from the genomic library to the gene pool.

        aGC's are converted to gGC's.
        Higher layer fields are updated.
        Nodes & edges are added the the GP graph.

        NOTE: This *MUST* be the only function pulling GC's into the GP from the GL.

        Args
        ----
        signatures (iterable(bytes[32])): Signatures to pull from the genomic library.
        """
        if _LOG_DEBUG: _logger.debug(f'Recursively pulling {signatures} into Gene Pool.')
        gcs = self._gl.recursive_select(_SIGNATURE_SQL, literals={'matches': signatures})
        self._gl.hl_copy(gcs)
        ggcs = {gc['ref']: gc for gc in map(partial(gGC, modified=True), gcs)}
        self.pool.update(ggcs)
        self.add_nodes(ggcs)

    def push_to_gp(self):
        """Insert or update locally modified gGC's into the persistent gene_pool.

        NOTE: This *MUST* be the only function pushing GC's to the persistent GP.
        """
        # TODO: This can be optimised to further minimise the amount of data munging of unmodified values.
        # TODO: Check for dodgy values that are not just bad logic e.g. overflows
        modified_gcs = [gc for gc in filter(_MODIFIED_FUNC, self.pool.values())]
        # FIXME: Use excluded columns depending on new or modified.
        for updated_gc in self._pool.upsert(modified_gcs, self._update_str, {}, _UPDATE_RETURNING_COLS):
            gc = self.pool[updated_gc['ref']]
            gc.update(updated_gc)
            for col in HIGHER_LAYER_COLS: # FIXME: Wrong definition - should be GP higher layer cols & use hl_copy().
                gc[col] = gc[col[1:]]
            gc['modified']= False

    def delete(self, refs):
        """Delete GCs from the gene pool.

        If ref is not in the pool it is ignored.
        If ref is aliased by a signature both the signature & ref entries are deleted.

        Args
        ----
        refs(iterable(int)): GC 'ref' values to delete.
        """
        #FIXME: This needs to recursively delete now unreferenced GC's
        refs = tuple(filter(lambda x: x in self.pool, refs))
        for ref in refs:
            gc = self.pool[ref]
            if gc.get('signature', None) is not None:
                del self.pool[gc['signature']]
            if gc.get('exec', None) is not None:
                remove_callable(gc)
            del self.pool[ref]
        if refs:
            # FIXME: How does this work if another worker is using the GC?
            self._pool.delete('{ref} in {ref_tuple}', {'ref_tuple': refs})

    def individuals(self, identifier=0):
        """Return a generator of individual gGCs for a population.

        The population is identified by identifier.

        Args
        ----
        identifier (str or int): Either the population name or index.

        Returns
        -------
        generator(gGC) of individuals in the population.
        """
        if isinstance(identifier, str):
            identifier = [p for p in self._population_data if p['name'] == identifier][0]['idx']
        return (gc for gc in filter(lambda x: x['population'] == identifier, self.pool.values()))

    def cull_target(self, pop_idx, ngen):
        """Reduce the target population to the population size.

        Individuals added to the GP by other workers are pulled in for the cull.
        This keeps workers populations aligned.

        Algorithm (in order):
            1. Keep the top _FITTEST_FRACTION fittest
            2. Keep the top _EVO_FRACTION most evolvable
            3. Keep the top _DIVERSITY_FRACTION most diverse* fitness scores
            4. Keep the top _UNIQUE_FRACTION most diverse structures
            5. Randomly select from the rest weighted by fitness.

        *_FRACTION is of the population size (truncated if a fraction of an individual)
        1-4 may be overlapping sets.
        """
        # TODO: Implement elitist functions
        population = list(self.individuals(pop_idx))
        exclusions = [i['ref'] for i in population]
        competition = tuple(self._gl.select(_COMPETITION_SQL, {'pop_idx': pop_idx, 'exclusions': exclusions}, ('fitness', 'ref')))
        population.extend(competition)
        population.extend((offspring[0] for offspring in ngen))
        if _LOG_DEBUG:
            _logger.debug(f"{len(competition)} remote individuals identified as a competition.")

        # Do the math
        weights = array([1.0 - i['fitness'][0] for i in population], float32)
        weights /= sum(weights)
        marked = len(population) - self.config['size']
        victims = [individual['ref'] for individual in choice(population, (marked,), False, weights)]

        # Delete any victims that are in the gene pool (NB: those from other workers are not so are ignored)
        self.delete(victims)

        # Pull any non-victims that were defined by other workers into the local GP cache.
        competition_winners = [i['ref'] for i in population if i['ref'] not in victims and 'f_count' not in i]
        if _LOG_DEBUG:
            _logger.debug(f"{len(competition_winners)} remote individuals pulled into local population.")
        self.pull_from_gl(competition_winners)

        # Add any new generation winners that were not cull victims to the GP
        ngen_winners = [i for i in ngen if i[0]['ref'] not in victims]
        self.add_to_gp_cache(ngen_winners)

        # Remove the callable of any ngen that did not make it.
        # Be aware that some sub-GC's may be in the local gene pool cache for other reasons.
        for gc_list in filter(lambda x: x[0]['ref'] in victims, ngen):
            for gc in filter(lambda x: x['ref'] not in self.pool, gc_list):
                remove_callable(gc)

    def cull_physical(self):
        """Remove any pGC's that are not active in any layer.

        A pGC that is not active in a layer has an False f_valid.
        """
        purge = False
        for layer, _ in filter(lambda x: x[1], enumerate(self.layer_evolutions)):
            cull_physical(self.pool.values(), layer)
            purge = True

        if purge:
            self.delete((i['ref'] for i in self.individuals(0) if not any(i['f_valid'])))

    def epoch(self, pop_idx):
        """Evolve the population and assess its fitness.

        Evolutionary steps are:
            1. Select a pGC to operate on the individual
            2. Evolve the individual into an offspring
            3. Assess the childs fitness
            4. Update the individuals parameters
            5. Update the parameters of the pGC
            6. Update the parameters of the pGC's pGC recursively

        Args
        ----
        pop_idx (int): The index of the target population to evolve
        fitness (callable): The fitness function of the target population.

        Returns
        -------
        ngen ([[eGC]]): The offsprings (evolved version) of the population & sub-GC's.
        """
        # pGC fitness updates happen en-masse after each layer has selected mutations.
        start = time()
        selection = [(individual, select_pGC(self, individual, 0)) for individual in self.individuals(pop_idx)]
        ngen = []
        for individual, pgc in selection:
            offspring = pgc.exec((individual,))
            xGC_inherit(offspring[0], individual, pgc)
            ngen.append(offspring)
            new_fitness = offspring[0]['fitness'][0] = fitness(offspring[0])
            delta_fitness = new_fitness - individual['fitness'][0]
            xGC_evolvability(individual, delta_fitness, 0)
            pGC_fitness(self, pgc, delta_fitness)
        self.cull_target(pop_idx, ngen)
        self.cull_physical()
        self.push_to_gp()
        self.metrics(time()-start)

    def metrics(self):
        """Calculate and record metrics."""
        self.layer_metrics()
        self.pgc_metrics()
        self.gp_metrics()
        self.layer_evolutions = [0] * self.max_depth

    def layer_metrics(self, duration):
        """Per layer metrics."""
        for evolutions, layer in filter(lambda x: x[0], enumerate(self.layer_evolutions)):
            fitness = []
            evolvability = []
            generation = []
            gc_count = []
            c_count = []
            filter_func = lambda x: x['f_valid']
            for individual in filter(filter_func, self.pool.values()):
                 fitness.append(individual['fitness'][layer])
                 evolvability.append(individual['evolvability'][layer])
                 generation.append(individual['generation'][layer])
                 gc_count.append(individual['gc_count'])
                 c_count.append(individual['c_count'])
            self._metrics['layer_metrics'].insert({
                'layer': layer,
                'f_max': max(fitness),
                'f_mean': mean(fitness),
                'f_min': min(fitness),
                'e_max': max(evolvability),
                'e_mean': mean(evolvability),
                'e_min': min(evolvability),
                'g_max': max(generation),
                'g_mean': mean(generation),
                'g_min': min(generation),
                'gcc_max': max(gc_count),
                'gcc_mean': mean(gc_count),
                'gcc_min': min(gc_count),
                'cc_max': max(c_count),
                'cc_mean': mean(c_count),
                'cc_min': min(c_count),
                'eps': evolutions / duration,
                'tag': 0,
                'worker_id': 0})

    def pgc_metrics(self, duration):
        """Per pGC layer metrics."""
        for evolutions, layer in filter(lambda x: x[0], enumerate(self.layer_evolutions)):
            ne = []
            filter_func = lambda x: x['f_valid']
            for individual in filter(filter_func, self.pool.values()):
                 fitness.append(individual['fitness'][layer])
                 evolvability.append(individual['evolvability'][layer])
                 generation.append(individual['generation'][layer])
                 gc_count.append(individual['gc_count'])
                 c_count.append(individual['c_count'])
            self._metrics['layer_metrics'].insert({
                'layer': layer,
                'f_max': max(fitness),
                'f_mean': mean(fitness),
                'f_min': min(fitness),
                'e_max': max(evolvability),
                'e_mean': mean(evolvability),
                'e_min': min(evolvability),
                'g_max': max(generation),
                'g_mean': mean(generation),
                'g_min': min(generation),
                'gcc_max': max(gc_count),
                'gcc_mean': mean(gc_count),
                'gcc_min': min(gc_count),
                'cc_max': max(c_count),
                'cc_mean': mean(c_count),
                'cc_min': min(c_count),
                'eps': evolutions / duration,
                'tag': 0,
                'worker_id': 0})
