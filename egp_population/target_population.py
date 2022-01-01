"""Target Population Management"""

from egp_physics.ep_type import vtype
from numpy.random import choice
from numpy import array, float32, sum
from .population import population


class target_population(population):

    # TODO: Default gene_pool should be localhost
    def __init__(self, fitness, diversity, inputs, outputs, gp, size=100, vt=vtype.OBJECT):
        """Create a target (non-physical GC's) population of individuals in a gene pool."""
        super.__init__(fitness, diversity, inputs, outputs, gp, size, vt)

    def cull(self, eppop):
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
        weights = array([1.0 - i['fitness'][0] for i in eppop], float32)
        weights /= sum(weights)
        marked = len(eppop) - self.config['size']
        return [individual['ref'] for individual in choice(eppop, (marked,), False, weights)] if marked else []
