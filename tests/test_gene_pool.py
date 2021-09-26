"""Test the Gene Pool."""


from egp_population.gene_pool import default_config as gp_default_config
from egp_population.gene_pool import gene_pool
from egp_genomic_library import default_config as gl_default_config
from egp_genomic_library import genomic_library


# Create a genomic library
_GL_CONFIG = gl_default_config()
_GL_CONFIG['database']['dbname'] = 'test_db'
_GL_CONFIG['delete_table'] = True
_GL = genomic_library(_GL_CONFIG)


# Gene pool config
_GP_CONFIG = gp_default_config()
_GP_CONFIG['database']['dbname'] = 'test_db'
_GP_CONFIG['delete_table'] = True


def test_default_instanciation():
    """Simple instanciation."""
    gene_pool(_GL, _GP_CONFIG)


def test_initialisation_1():
    """Initialise the gene pool with a codon interface definition.

    Using a codon interface definition guarantees a match.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.initialize((0, 0), (0,))
    assert len(gp._pool) == 1


def test_initialisation_3():
    """Initialise the gene pool with a non-codon interface definition.

    Using a codon interface definition guarantees a match.
    """
    gp = gene_pool(_GL, _GP_CONFIG)
    gp.initialize((0, 0, 0.0, 0.0), (0.0,), num=5)
    assert len(gp._pool) == 5
