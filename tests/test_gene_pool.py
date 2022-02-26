"""Test the Gene Pool."""

from logging import NullHandler, getLogger
from egp_population.gene_pool import default_config as gp_default_config
from egp_population.gene_pool import gene_pool
from egp_genomic_library import default_config as gl_default_config
from egp_genomic_library import genomic_library
from egp_physics.ep_type import vtype
from numpy.random import randint
from numpy import int32, float32, where, arange, array, clip, isfinite

# Logging
_logger = getLogger(__name__)
_logger.addHandler(NullHandler())

# Create a genomic library
_GL_CONFIG = gl_default_config()
_GL_CONFIG['database']['dbname'] = 'test_db'
_GL_CONFIG['delete_table'] = True
_GL = genomic_library(_GL_CONFIG)


# Gene pool config
_GP_CONFIG = gp_default_config()
_GP_CONFIG['gp']['database']['dbname'] = 'test_db'
for table in _GP_CONFIG:
    _GP_CONFIG[table]['delete_table'] = True

# Divide! Characterization function
_X1_TEST = randint(-100,100,(100,),dtype=int32).astype(float32)
_X2_TEST = randint(-100,100,(100,),dtype=int32).astype(float32)
_X2_TEST = where(_X2_TEST == 0, arange(100, dtype=float32), _X2_TEST)
_Y_TEST = _X1_TEST /_X2_TEST
def characterize(gc):
    """"Characterize gc for Divide!

    Fitness is 1 - the normalised clipped mean-squared-error of 100 random division examples.
    The squared error is clipped to 10 and normalised to 1.

    Survivability is fitness in this test.

    NOTE: GGC may be modified in a limited way

    Args
    ----
    gc (gGC): The GC to be characterized.

    Returns
    -------
    tuple(float, float): Fitness, Survivability both in the range 0.0 <= x <= 1.0
    """
    results = []
    for x1, x2, y in zip(_X1_TEST, _X2_TEST, _Y_TEST):
        result = gc['exec']((x1, x2))
        if result is None:
            result = [100000.0]
        results.append(result[0])
    y_pred = array(results, dtype=float32)
    y_pred = where(isfinite(y_pred), y_pred, 100000.0)
    _logger.debug(f'GC {gc["ref"]} Predicted: {y_pred}')
    clipped = clip((_Y_TEST - y_pred) ** 2, -10.0, 10.0)
    _logger.debug(f'GC {gc["ref"]} Clipped: {clipped}')
    mse = clipped.mean() / 10
    _logger.debug(f'GC {gc["ref"]} MSE: {mse}')
    fitness = (1.0 - mse).sum()
    _logger.debug(f'GC {gc["ref"]} fitness = {fitness}, survivability = {fitness}')
    return fitness, fitness

def recharacterize(gcs):
    """Re-characterize gc for Divide!

    The 'survivability'[0] value of each GC is modified.

    Args
    ----
    gcs (tuple(gGC)): The GCs to be recharacterized.
    """
    # In this case no modification is made.
    pass

# A population config
_P_CONFIG = {
    'size': 100,
    'name': 'Divide!',
    'inputs': ('float', 'float'),
    'outputs': ('float',),
    'characterize': characterize,
    'recharacterize': recharacterize,
    'vt': vtype.EP_TYPE_STR,
    'description': 'Input values are x1, x2. Desired return value, y, is x1/x2.'
}

def test_default_instanciation():
    """Simple instanciation."""
    gene_pool(_GL, _GP_CONFIG)


def test_initialisation():
    """Initialise the gene pool with a codon interface definition.

    Using a codon interface definition guarantees a match.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.create_population(_P_CONFIG)


def test_evolve_simple():
    """Initialise the gene pool with a non-codon interface definition.

    Forces new GC's to be generated.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.create_population(_P_CONFIG)
    gp.evolve(num_sub_processes=1)
