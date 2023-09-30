"""Population Types"""
from datetime import datetime
from typing import Callable, NotRequired, TypedDict, Required
from uuid import UUID

from numpy import int64, single, bool_
from numpy.typing import NDArray

from .population import population


Survivability = NDArray[single]
Active = NDArray[bool_]
Reference = NDArray[int64]
FitnessFunction = Callable[[Callable], single]
SurvivabilityFunction = Callable[
    [population],
    tuple[Survivability, Active] | tuple[Reference, Survivability, Active],
]


class ProblemConfig(TypedDict):
    """Type definition."""

    name: NotRequired[str | None]
    description: NotRequired[str | None]
    inputs: Required[list[str]]
    outputs: Required[list[str]]
    creator: NotRequired[UUID | None]


class PopulationConfig(TypedDict):
    """Type definition."""

    uid: NotRequired[int]
    egp_problem: NotRequired[str]
    worker_id: NotRequired[UUID | str]
    size: NotRequired[int]
    name: Required[str]
    inputs: NotRequired[list[str]]
    outputs: NotRequired[list[str]]
    ordered_interface_hash: NotRequired[int]
    unordered_interface_hash: NotRequired[int]
    description: NotRequired[str | None]
    meta_data: NotRequired[str | None]
    created: NotRequired[datetime | str]
    updated: NotRequired[datetime | str]
    survivability: NotRequired[str]
    survivability_function: NotRequired[SurvivabilityFunction]
    fitness_function: NotRequired[FitnessFunction]


class PopulationConfigNorm(TypedDict):
    """Type definition."""

    uid: int
    egp_problem: str
    worker_id: UUID
    size: int
    inputs: list[str]
    outputs: list[str]
    ordered_interface_hash: int
    unordered_interface_hash: int
    name: str
    description: str | None
    meta_data: str | None
    created: datetime
    updated: datetime
    survivability_function: SurvivabilityFunction
    fitness_function: FitnessFunction


class PopulationsConfig(TypedDict):
    """Type definition."""

    worker_id: NotRequired[UUID | str]
    configs: NotRequired[list[PopulationConfig]]


class PopulationsConfigNorm(TypedDict):
    """Type definition."""

    worker_id: UUID
    configs: list[PopulationConfigNorm]
