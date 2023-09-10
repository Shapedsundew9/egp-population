"""Erasmus GP populations module."""
from functools import lru_cache
from logging import DEBUG, Logger, NullHandler, getLogger
from random import getrandbits
from typing import Any, Generator, Self

from egp_types.xGC import xGC
from numpy import ndarray, array

_logger: Logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG: bool = _logger.isEnabledFor(DEBUG)


class population:
    """The population class provides an interface to the population in the
    Gene Pool typically for use in the survivor selection process."""

    def __init__(
        self, xgcs: list[xGC] | tuple[xGC, ...] | Generator[xGC, None, None]
    ) -> None:
        """Create a new population instance.

        NOTE: to prevent copying and hashing of the xGC's the xGC's are stored in a tuple.
        """
        self._xgcs: tuple[xGC, ...] | list[xGC] = (
            tuple(xgcs) if not isinstance(xgcs, (tuple, list)) else xgcs
        )
        self._hash: int = getrandbits(64)

    def __eq__(self, other: Any) -> bool:
        """Required for LRU cache to work for each instance."""
        if isinstance(other, population):
            return self._hash == other._hash
        return False

    def __hash__(self) -> int:
        """Required for LRU cache to work for each instance."""
        return self._hash

    def __len__(self) -> int:
        """The number of xGCs in the population."""
        return len(self._xgcs)

    def active(self) -> Self:
        """Return the population of active xGCs."""
        return population(xgc for xgc in self._xgcs if xgc["active"])

    def modified(self) -> None:
        """If the xGC list changes or the undelrying data changes then change the hash.

        This is necessary because it would be expensive to check the xGC's for changes.
        """
        self._hash = getrandbits(64)

    @lru_cache(maxsize=64)
    def __getitem__(self, key: str) -> list | ndarray:
        """Return the list of values for the given key."""
        # NOTE: The LRU cache is class level so it will be shared between all instances.
        return array([xgc[key] for xgc in self._xgcs])
