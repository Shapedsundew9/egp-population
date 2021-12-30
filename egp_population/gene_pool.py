"""Gene pool management for Erasmus GP."""

from copy import deepcopy
from functools import partial
from os.path import dirname, join
from logging import DEBUG, INFO, WARN, ERROR, FATAL, NullHandler, getLogger
from json import load, loads, dumps
from random import random
from numpy.core.fromnumeric import mean
from numpy.random import choice
from numpy import array, float32, sum, zeros
from egp_genomic_library.genetic_material_store import genetic_material_store
from egp_genomic_library.genomic_library import (compress_json, decompress_json, UPDATE_STR, UPDATE_RETURNING_COLS,
    sha256_to_str, sql_functions, HIGHER_LAYER_COLS)
from egp_physics.ep_type import asstr, vtype
from egp_physics.gc_type import eGC, interface_definition, gGC, interface_hash
from egp_physics.physics import stablise, xGC_evolvability, pGC_fitness, select_pGC, xGC_inherit, random_reference
from egp_physics.gc_graph import gc_graph
from time import time, sleep


from pypgtable import table
from .gpm import create_callable, remove_callable


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
_GENE_POOL_MODULE = 'gpm'
_GPP_COLUMNS = ('idx', 'size', 'name', 'description', 'meta_data')
_POPULATION_IN_DEFINITION = -1


# The physical GC population
_PHYSICAL_ENVIRONMENT = {
    'fitness': fitness,
    'diversity': diversity,
    'size': 0,
    'inputs': ('gGC',),
    'outputs': ('gGC',),
    'vt': vtype.INSTANCE_STR
}

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
with open(join(dirname(__file__), "formats/gpp_table_format.json"), "r") as file_ptr:
    _GPP_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/gp_metrics_format.json"), "r") as file_ptr:
    _GP_METRICS_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/layer_metrics_format.json"), "r") as file_ptr:
    _LAYER_METRICS_TABLE_SCHEMA = load(file_ptr)
with open(join(dirname(__file__), "formats/pgc_metrics_format.json"), "r") as file_ptr:
    _PGC_METRICS_TABLE_SCHEMA = load(file_ptr)


# GC queries
_SIGNATURE_SQL = 'WHERE ({signature} = ANY({matches}))'
_INITIAL_GC_SQL = ('WHERE {input_types} = {itypes} AND {inputs}={iidx} '
                      'AND {output_types} = {otypes} AND {outputs}={oidx} '
                      'AND NOT({exclude_column} = ANY({exclusions})) ORDER BY RANDOM() LIMIT {limit}')


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
    mapping both 'ref' and 'signature' keys to gGC's.
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

        # pool is a dictionary of ref:gGC and signature:gGC
        # All gGC's in pool must be valid & any modified are sync'd to the gene pool table
        # in the database at the end of the target epoch.
        self.pool = {}
        self._gl = genomic_library
        self._pool = table(configs['gp'])

        # UID for this worker
        self.worker_id = random_reference()

        # FIXME: Validators for all tables.

        # A dictionary of metric tables. Metric tables are undated at the end of the
        # taregt epoch.
        self._metrics = {c: table(configs[c]) for c in configs.keys() if 'metrics' in c}

        # ?
        self._env_config = None
        self._populations_table = self._create_populations_table(configs['gp'])

        # Modify the update strings to use the right table for the gene pool.
        self._update_str = UPDATE_STR.replace('__table__', configs['gp']['table'])

        #?
        self._interface = None
        self.encode_value = self._pool.encode_value
        self.select = self._pool.select

        # The number of layers in the gene pool
        # Layer 0 is the tsrget layer
        # Layers 1 to N are pGC layers
        self.num_layers = 0

        # Used to track the number of updates to individuals in a layer (?)
        self.layer_evolutions = [0]

        # If this instance created the gene pool then it is responsible for configuring
        # setting up database functions and initial population.
        if self._pool.raw.creator:
            self._pool.raw.arbitrary_sql(sql_functions(), read=False)


    def _create_populations_table(self, config):
        """Create a table to store population information.

        The populations table has all the same DB details as the gene
        pool.

        Args
        ----
        config (gene_pool config): The config of the gene pool

        Returns
        -------
        pypgtable.table
        """
        _config = {'database': config['database']}
        _config['table'] = config['table'] + '_populations'
        _config['create_table'] = True
        _config['schema'] = _GPP_TABLE_SCHEMA
        tbl = table(_config)
        tbl.register_conversion('inputs', dumps, loads)
        tbl.register_conversion('outputs', dumps, loads)
        return tbl

    def populations(self):
        """Return the definition of all populations."""
        return self._populations_table.select()

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
        1.0s (on average) until size changes from _POPULATION_IN_DEFINITION.

        Args
        ----
        config (dict): (required)
            'size':(int): Defines the population size.
            'name':(str): An unique short (<= 64 charaters) string giving a name.
            'inputs':iter(vtype): Target individual input definition using vt type.
            'outputs':iter(vtype): Target individual output definition using vt type.
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
        data['vt'] = vtype.EP_TYPE_STR
        data = next(self._populations_table.insert((data,), '*'))

        # If this worker created the entry then populate the gene pool.
        if data['worker_id'] == self.worker_id:
            self.populate(data['inputs'], data['outputs'], data['idx'], num=config['size'], vt=data['vt'])
            data['size'] = config['size']
            self._populations_table.update("{size} = {s}", "{idx} = {i}", {'s': data['size'], 'i': data['idx']})
        else:
            while data['size'] == _POPULATION_IN_DEFINITION:
                _logger.info(f"Waiting for population {data['name']} to be created by worker {data['worker_id']}.")
                sleep(0.5 + random())
                data = next(self._populations_table.select('{name} = {n}', {'n': data['name']}))
        return data

    def populate(self, inputs, outputs, pop_idx, exclusions=[], num=1, vt=vtype.OBJECT):
        """Fetch or create num valid GC's with the specified inputs & outputs.

        The construct a population with inputs and outputs as defined by inputs, outputs & vt.

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
        pop_idx (int): The population index to use for the new individuals.
        num (int): The number of GC's to create (or delete if -ve)
        vt (vtype): See vtype definition.

        Returns
        -------
        population index (int) or None.
        """
        # Check the population index exists
        literals = {'pop_idx': pop_idx}
        result = next(self._populations.select('WHERE {idx} = {pop_idx}', literals))
        if not len(result):
            raise ValueError(f"Population index {pop_idx} was not found.")
        count = len(tuple(self._pool.raw.select('WHERE {population} = {pop_idx}', literals, columns=('population',))))
        num = result['size'] - count
        verb, tf = ('Adding', 'to') if num > 0 else ('Deleting', 'from')
        _logger.info(f"{verb} {abs(num)} GC's {tf} population index: {pop_idx}, name: {result['name']}")

        # Create search criteria for the genomic library
        xputs = {
            'exclude_column': 'signature',
            'exclusions': exclusions,
            'limit': num
        }
        input_eps, xputs['itypes'], xputs['iidx'] = interface_definition(inputs, vt)
        output_eps, xputs['otypes'], xputs['oidx'] = interface_definition(outputs, vt)
        self._interface = interface_hash(input_eps, output_eps)

        # Adding to the gene pool
        if num > 0:
            # Find the GC's that match and then recursively pull them from the genomic library
            matches = list(m[0] for m in self._gl.library.raw.select(_INITIAL_GC_SQL, literals=xputs, columns=('signature',)))
            _logger.info(f'Found {len(matches)} GCs matching input and output criteria.')
            if _LOG_INFO and matches:
                _logger.debug(f'GC signatures: {[sha256_to_str(m) for m in matches]}.')
            self.pull(matches)
            for gc in filter(lambda x: x['signature'] in matches, self.pool.values()):
                gc['modified'] = True
                gc['population'] = pop_idx

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
            self.add(ggcs)

        # Removing from the gene pool
        elif num < 0:
            xputs['limit'] = -num
            matches = list(m[0] for m in self._pool.raw.select(_INITIAL_GC_SQL, literals=xputs, columns=('signature',)))
            self.delete(matches)

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
            if _LOG_DEBUG: _logger.debug(f'Added {ggc["ref"]} to local cached gene pool.')
            gca = ggc.get('gca', None)
            if gca is not None and ggc['gca_ref'] not in self.pool:
                children.append(gca)
                if _LOG_DEBUG: _logger.debug(f'Appended {gca} to list to pull from the Genomic Library.')
            gcb = ggc.get('gcb', None)
            if gcb is not None and ggc['gcb_ref'] not in self.pool:
                children.append(gcb)
                if _LOG_DEBUG: _logger.debug(f'Appended {gcb} to list to pull from the Genomic Library.')
        self.pull(children)
        self.add_nodes((ggc for ggc in filter(_MODIFIED_FUNC, self.pool.values())))
        #FIXME: Must correctly update higher layer fields & not wipe out this layer.
        #TODO: Not sure if this is the right place to push. Might need to update
        # locally on every change to maintain consistency but do not what to do a DB push everytime.
        self.push()

    def pull(self, signatures):
        """Pull aGCs and all sub-GC's recursively from the genomic library to the gene pool.

        aGC's are converted to gGC's.

        Args
        ----
        signatures (iterable(bytes[32])): Signatures to pull from the genomic library.
        """
        if _LOG_DEBUG: _logger.debug(f'Recursively pulling {signatures} into Gene Pool.')
        gcs = self._gl.recursive_select(_SIGNATURE_SQL, literals={'matches': signatures})
        self.pool.update({gc['ref']: gc for gc in map(partial(gGC, modified=True), gcs)})

    def push(self):
        """Insert or update into locally modified gGC's into the persistent gene_pool."""
        # TODO: This can be optimised to further minimise the amount of data munging of unmodified values.
        modified_gcs = [gc for gc in filter(_MODIFIED_FUNC, self.pool.values())]
        for updated_gc in self._pool.upsert(modified_gcs, self._update_str, {}, _UPDATE_RETURNING_COLS):
            gc = self.pool[updated_gc['ref']]
            gc.update(updated_gc)
            for col in HIGHER_LAYER_COLS:
                gc[col] = gc[col[1:]]
            gc['modified']= False
        self.define(modified_gcs)

    def define(self, gcs):
        """Define the executable object for each gc in gcs."""
        for gc in gcs:
            if gc.get('exec', None) is None:
                create_callable(gc)

    def undefine(self, gcs):
        """Delete the executable object for each gc in gcs."""
        for gc in gcs:
            if gc.get('exec', None) is not None:
                remove_callable(gc)

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
            self._pool.delete('{ref} in {ref_tuple}', {'ref_tuple': refs})

    def individuals(self, pop_idx=0):
        """Return a generator of individual gGCs."""
        return (gc for gc in filter(lambda x: x['population'] == pop_idx, self.pool.values()))

    def cull_target(self, pop_idx=0):
        """Reduce the target population to the population size.

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
        population = tuple(self.individuals(pop_idx))
        weights = array([1.0 - i['fitness'][0] for i in population], float32)
        weights /= sum(weights)
        marked = len(population) - self.config['size']
        self.delete([individual['ref'] for individual in choice(population, (marked,), False, weights)] if marked else [])

    def evolve_target(self, pop_idx, fitness):
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
            self.define(offspring)
            new_fitness = offspring[0]['fitness'][0] = fitness(offspring[0])
            delta_fitness = new_fitness - individual['fitness'][0]
            xGC_evolvability(individual, delta_fitness, 0)
            pGC_fitness(self, pgc, delta_fitness)
        self.layer_evolutions[0] = len(selection)
        self.cull_target()
        self.push()
        # Pull Individuals that have been created by other workers
        # Cull again
        self.metrics(time()-start)
        return ngen

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
            filter_func = lambda x:len(x['f_count']) > layer and x['f_count']
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
            filter_func = lambda x:len(x['f_count']) > layer and x['f_count']
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
