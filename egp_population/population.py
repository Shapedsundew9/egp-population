"""Population Management"""

from .gene_pool import gene_pool
from egp_physics.ep_type import vtype


class population():


    def __init__(self, fitness, inputs, outputs, gp: gene_pool, num=100, vt=vtype.OBJECT):
        """Create a population of individuals in a gene pool.

        The gene pool may be empty
        """
        self.gp = gp
        self.gp.initialize(inputs, outputs, num=num, vt=vt)
        self.fitness = fitness
