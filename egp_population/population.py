"""Population Management"""

from egp_physics.ep_type import vtype
from numpy.random import choice
from numpy import array, float32, sum


#Cull algorithm parameters
_FITTEST_FRACTION = 0.03
_EVO_FRACTION = 0.05
_DIVERSE_FRACTION = 0.05
_UNIQUE_FRACTION = 0.05
assert _FITTEST_FRACTION + _EVO_FRACTION + _DIVERSE_FRACTION + _UNIQUE_FRACTION < 1.0


# No size limit population contraints
MINIMUM_FITNESS = 0.01
MINIMUM_POPULATION = 100


class population():

    # TODO: Default gene_pool should be localhost
    def __init__(self, config):
        """Create a population of individuals in a gene pool.

        Args
        ----
        config (dict):
            'fitness' (callable): A function that takes an individual GC callable ('exec') and
                returns a tuple of (fitness score, results) where
                    0.0 <= fitness score <= 1.0
                    results is a tuple of individual test case results (objects)
            'diversity' (callable): A function that takes a list/tuple of results (see fitness) and
                returns a list/tuple of indices of results that make up the smallest most diverse set.
            'inputs' (iterable(object or type)): Input objects of vt type.
            'outputs' (iteratble(object or type)): Output objects of vt type.
            'pop_idx' (int): The population index to use for the new individuals.
            'size' (int): The size of the population. If 0 growth of the population is constrained to
                MIN_POPULATION individuals or individuals with >= MINIMUM_FITNESS long term fitness
                whichever is greater.
            'depth' (int): The layer of the population in the environment. Target populations must have a
                depth of 0. depth must be >= 0.
            'vt' (vtype): See vtype definition.
        """
        self.gp = config['gp']
        self.config = self.gp.upsert_population(config)
        self.depth = config['depth']
        self.fitness = config['fitness']
        self.diversity = config['diversity']

    def epoch(self):
        """Age the population by 1 generation.

        Within the gene pool there are physical GC's that modifiy the population individuals.
        The physical GC's take one individual and return a modified/different individual.
        """
        # new_population is a generator of iterables of eGCs with the first being the top level GC
        # of an individuals child
        new_ggcs = []
        new_population = []
        for individual in self.gp.individuals(self.config['idx']):
            child = self.gp.evolve(individual, self)
            new_ggcs.append(child)
            new_population.append(child[0])
        new_population.extend(self.gc.individuals(self.config['idx']))

        # Trim the new population to size
        marked = self.cull(new_population)
        self.gp.delete(marked)
        ggcs = []
        for ggc in filter(lambda x: x[0]['ref'] not in marked, new_ggcs):
            ggcs.extend(ggc)
        self.gc.add(ggcs)

    def cull(self, eppop, depth=0):
        """Reduce eppop to the population size.

        If the population size limit is 0 then the population is constrained to
        MIN_POPULATION individuals or individuals with >= MINIMUM_FITNESS long term fitness.

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
        weights = array([1.0 - i['fitness'][depth] for i in eppop], float32)
        weights /= sum(weights)
        marked = len(eppop) - self.config['size']
        return [individual['ref'] for individual in choice(eppop, (marked,), False, weights)] if marked else []