"""Test the Gene Pool."""
from warnings import filterwarnings

filterwarnings("ignore", category=DeprecationWarning)

from logging import NullHandler, getLogger
from egp_population.gene_pool import default_config as gp_default_config
from egp_population.gene_pool import gene_pool
from egp_genomic_library import default_config as gl_default_config
from egp_genomic_library import genomic_library
from egp_physics.ep_type import vtype
from numpy.random import randint, normal
from numpy import int32, float32, where, arange, array, clip, isfinite, isnan
from scipy.stats import linregress

# Logging
_logger = getLogger(__name__)
_logger.addHandler(NullHandler())

# Create a genomic library
_GL_CONFIG = gl_default_config()
_GL_CONFIG["database"]["dbname"] = "test_db"
_GL_CONFIG["delete_table"] = True
_GL = genomic_library(_GL_CONFIG)


# Gene pool config
_GP_CONFIG = gp_default_config()
for table in _GP_CONFIG:
    _GP_CONFIG[table]["delete_table"] = True
    _GP_CONFIG[table]["database"]["dbname"] = "test_db"


# Divide! Characterization function
_X1_TEST = randint(-100, 100, (100,), dtype=int32).astype(float32)
_X2_TEST = randint(-100, 100, (100,), dtype=int32).astype(float32)
_X2_TEST = where(_X2_TEST == 0, arange(100, dtype=float32), _X2_TEST)
_Y_TEST = _X1_TEST / _X2_TEST


def characterize_divide(gc):
    """ "Characterize gc for Divide!

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
        result = gc["exec"]((x1, x2))
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


def recharacterize_divide(gcs):
    """Re-characterize gc for Divide!

    The 'survivability'[0] value of each GC is modified.

    Args
    ----
    gcs (tuple(gGC)): The GCs to be recharacterized.
    """
    # In this case no modification is made.
    pass


# A population config
_DIV_CONFIG = {
    "size": 100,
    "name": "Divide!",
    "inputs": ("float", "float"),
    "outputs": ("float",),
    "characterize": characterize_divide,
    "recharacterize": recharacterize_divide,
    "vt": vtype.EP_TYPE_STR,
    "description": "Input values are x1, x2. Desired return value, y, is x1/x2.",
}

# Linear Fit Characterization function
_M = 17.42
_C = -4.2
_MSE_LIMIT = 2.0**12
_MIN_FITNESS = 1.0 / _MSE_LIMIT


def characterize_fit(gc):
    """ "Characterize gc for a linear function fit!

    Fitness is means squared error of 101 with unit separation in the range
    -50 to +50

    Survivability is fitness in this test.

    NOTE: GGC may be modified in a limited way

    Args
    ----
    gc (gGC): The GC to be characterized.

    Returns
    -------
    tuple(float, float): Fitness, Survivability both in the range 0.0 <= x <= 1.0 or None
    """
    results = [gc["exec"]((x)) for x in arange(-50, 51)]
    if None in results:
        fitness = 0.0
    else:
        mse = ((results - (arange(-50, 51) * _M + _C)) ** 2).mean()
        if not isfinite(mse):
            fitness = 0.0
        elif mse > _MSE_LIMIT:
            fitness = _MIN_FITNESS
        else:
            fitness = 1.0 - mse / _MSE_LIMIT
    _logger.debug(f'GC {gc["ref"]} fitness = {fitness}, survivability = {fitness}')
    return fitness, fitness


def recharacterize_fit(gcs):
    """Re-characterize gc for linear fit!

    The 'survivability'[0] value of each GC is modified.

    Args
    ----
    gcs (tuple(gGC)): The GCs to be recharacterized.
    """
    # In this case no modification is made.
    pass


# A population config
_FIT_CONFIG = {
    "size": 100,
    "name": "Linear fit",
    "inputs": ("float",),
    "outputs": ("float",),
    "characterize": characterize_fit,
    "recharacterize": recharacterize_fit,
    "vt": vtype.EP_TYPE_STR,
    "description": "Input values are x. Desired return value, y = _M * x + _C",
}


def test_default_instanciation():
    """Simple instanciation."""
    gene_pool(_GL, _GP_CONFIG)


def test_initialisation():
    """Initialise the gene pool with a codon interface definition.

    Using a codon interface definition guarantees a match.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.create_population(_DIV_CONFIG)


def test_evolve_simple_1():
    """Initialise the gene pool with a non-codon interface definition.

    Forces new GC's to be generated.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.create_population(_DIV_CONFIG)
    _logger.debug("Population created.")
    gp.evolve(num_sub_processes=0)


def test_evolve_simple_2():
    """Initialise the gene pool with a non-codon interface definition.

    Forces new GC's to be generated.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.create_population(_FIT_CONFIG)
    _logger.debug("Population created.")
    gp.evolve(num_sub_processes=0)
