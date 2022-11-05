"""Gene Pool Cache.

The gene pool cache is a space and time optimised store of GC's. It is designed to
be multi-process friendly.

Naively, the gene pool cache could be implemented as a dictionary with reference keys.
This would be fast but does not scale well. Python dictionaries use huge amounts
of memory and are updated in a spatially broad manner requiring subprocesses to maintain
an almost full copy even if most entries are only read.

The gene pool cache as implemented here maintains a dictionary like interface but takes
advantage of some GC structural design choices to efficiently store data in graph-tool
properties and fast index them with a trivial (at least as fast as a dictionary) lookup. 
It also takes advantage of GC usage patterns to cluster stable and volatile GC data which
makes efficient use of OS CoW behaviour in a multi-process environment as well as minimising
the volume of null data between GC variants by layering.
"""

from logging import DEBUG, NullHandler, getLogger
from typing import Self

# Logging
_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)


# Critical constants
# SP_UID_MSB: Sub-Process Unique IDentifier Most Significant Bit
# GPC_UID_MSB: Gene Pool Class Unique IDentifier Most Significant Bit
# GCI_UID_MSB: Genetic Code Index Unique IDentifier Most Significant Bit
_SP_UID_MSB = 63
_SP_UID_LSB = 32
_LAYER_MSB = 31
_LAYER_LSB = 30
_VERTEX_MSB = 29
_VERTEX_LSB = 0

# Masks for top level reference bit fields
_SPB_MASK = ((1 << (_SP_UID_MSB + 1)) - 1) ^ ((1 << _SP_UID_LSB) - 1)
_LAYER_MASK = ((1 << _SP_UID_LSB) - 1) ^ ((1 << _LAYER_LSB) - 1)
_VERTEX_MASK = (1 << _VERTEX_LSB) - 1

# Used in the Gene Pool to bound the sub-process UIDs
_SPB_UID_WIDTH = _SP_UID_MSB - _LAYER_MSB
SP_UID_LOWER_LIMIT = -(1 << (_SPB_UID_WIDTH - 1))
SP_UID_UPPER_LIMIT = (1 << (_SPB_UID_WIDTH - 1)) - 1

# Check the masks
_logger.debug(f'Sub-process UID mask: {_SPB_MASK:016X}')
_logger.debug(f'Layer mask: {_LAYER_MASK:016X}')
_logger.debug(f'Vertex mask: {_VERTEX_MASK:016X}')


# FIXME: This is temporary
# The Gene Pool Cache (GPC) has a dictionary like interface plus some additional
# access functions for better performance of common operations. An actual dictionary
# uses way too many resources but is easier to implement.
# In the short term (and may be in the long term as a validation reference)
# the GPC is implemented as a dictionary and the 'optimised' access
# functions emulated.

# NOTE the GPC is a global object. It is instanciated at the bottom of this file.

# The 'temporary' GPC
class gene_pool_cache(dict):
    pass


# The 'optimised' GPC
class _gpc():
    """Mapping to a GC stored in the gene pool cache."""

    __slots__ = '_gp', '_ref'

    def __init__(self, llgpc: _gpc | None):

