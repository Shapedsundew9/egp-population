""" The survivability function is used to determine which solutions continue to be evolved & which will be deleted.
    It is called after each epoch and has 2 ways to define the next epoch:
        1. Set 'active' to True for the individuals to continue to be evolved. Survivability is ignored.
        2. Set 'active' to False for all individuals and 'active_size' individuals with the highest survivability will be evolved.
    NOTE: In the case of ties in survivability crossing the 'active_size' limit the individuals with the same survivability
    will be chosen at random until the limit is reached.

    Erasmus does not delete individuals unless memory is exhausted (so the total population can grow to be very large).
    The intention is not to close off evolutionary paths that may ultimately lead to better solutions. When forced to
    delete individuals only ones with no existing descendants will be deleted. The ones with the lowest survivability first.
    In the case of ties individuals will be chosen at random.

    Trivially, and by default, survivability == fitness.
        0.0 <= survivability <= 1.0
    However, the power of survivability is that it allows a more complex definition of what constitutes the best solution without
    obfuscating the fitness function and can be evaluated in the context of the population rather than the just the individual.
    For example, solution diversity can be encouraged by increasing the survivabilty of individuals that are structurally
    different. In the early stages of evolution a sub-population may show increased fitness and start to dominate the active
    population but run into a deadend, successively failing to improve. If survivability were based only on an individuals
    performance a slightly less fit individual with a different structure and different potential may be deleted before it has a
    chance to be evolved and that time/energy/cost would be lost.

    Returns
    -------
    tuple[array[single], array[bool]] = (survivability, active) for the whole population
    tuple[array[int64], array[single], array[bool]] = (ref, survivability, active) for a subset of the population
"""
from typing import cast

from numpy import bool_, full, single
from numpy.typing import NDArray

from .population import population

from .egp_typing import Active, Reference, Survivability, SurvivabilityFunction


def elitist(populous: population) -> tuple[Survivability, Active] | tuple[Reference, Survivability, Active]:
    """Elitist survivability function. The fittest active individuals in the gene pool will be evolved."""
    survivability: NDArray[single] = cast(NDArray[single], populous["fitness"])
    active: NDArray[bool_] = full(survivability.shape, False, dtype=bool_)
    return (survivability, active)


SURVIVABILITY_FUNCTIONS: dict[str, SurvivabilityFunction] = {"elitist": elitist}
