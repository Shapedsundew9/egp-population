"""Test the Gene Pool."""

from logging import NullHandler, getLogger
from egp_population.gene_pool import default_config as gp_default_config
from egp_population.gene_pool import gene_pool
from egp_genomic_library import default_config as gl_default_config
from egp_genomic_library import genomic_library
from egp_physics.ep_type import vtype

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


# A population config
_P_CONFIG = {
    'size': 100,
    'name': 'Divide!',
    'inputs': ('float', 'float'),
    'outputs': ('float',),
    'characterize': lambda x: x,
    'vt': vtype.EP_TYPE_STR,
    'description': 'Input values are x, y. Desired return value is x/y.'
}

def test_default_instanciation():
    """Simple instanciation."""
    gene_pool(_GL, _GP_CONFIG)


def test_initialisation_1():
    """Initialise the gene pool with a codon interface definition.

    Using a codon interface definition guarantees a match.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.create_population(_P_CONFIG)
    assert len(gp._pool) == 1


def test_initialisation_3():
    """Initialise the gene pool with a non-codon interface definition.

    Forces new GC's to be generated.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    config = gp.upsert_population({'size': 100, 'inputs': (0, 0, 0.0, 0.0), 'outputs': (0.0,), 'vt': vtype.OBJECT})
    for individual in gp.individuals(config['idx']):
        _logger.debug(f"Executing individual ref {individual['ref']}:")
        try:
            retval = individual['exec']((1, 2, 1.0, 2.0))
        except (ZeroDivisionError, ValueError, IndexError) as e:
            _logger.warning(f"Individual ref {individual['ref']} threw exception '{str(e)}'")

        assert isinstance(retval, tuple)
        assert isinstance(retval[0], float)


if __name__ == "__main__":
    test_initialisation_3()