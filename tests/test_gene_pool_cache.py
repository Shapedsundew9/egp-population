"""Test the Gene Pool."""

from logging import NullHandler, getLogger
from egp_population.gene_pool_cache import gene_pool_cache

# Logging
_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
