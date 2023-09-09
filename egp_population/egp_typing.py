"""Population Types"""
from datetime import datetime
from typing import Callable, NotRequired, TypedDict
from uuid import UUID

from numpy import int64, single, bool_
from numpy.typing import NDArray

from .population import population


FitnessFunction = Callable[[Callable], single]
SurvivabilityFunction = Callable[[population], tuple[NDArray[single], NDArray[bool_]] | tuple[NDArray[int64], NDArray[single], NDArray[bool_]]]


class PopulationConfig(TypedDict):
    """Type definition."""
    uid: NotRequired[int]
    population_hash: NotRequired[bytes]
    git_repo: NotRequired[str]
    git_url: NotRequired[str]
    git_hash: NotRequired[str]
    verified: NotRequired[bool]
    worker_id: NotRequired[UUID]
    size: NotRequired[int]
    inputs: NotRequired[tuple[str, ...]]
    outputs: NotRequired[tuple[str, ...]]
    ordered_interface_hash: NotRequired[int]
    name: NotRequired[str]
    description: NotRequired[str | None]
    meta_data: NotRequired[str | None]
    created: NotRequired[datetime]
    updated: NotRequired[datetime]
    preload_function: NotRequired[Callable[[], None]]
    fitness_function: NotRequired[FitnessFunction]
    survivability_function: NotRequired[SurvivabilityFunction]
    preload_import: NotRequired[str]
    fitness_import: NotRequired[str]
    survivability_import: NotRequired[str]


class PopulationConfigNorm(TypedDict):
    """Type definition."""
    uid: int
    population_hash: bytes
    git_repo: str
    git_url: str
    git_hash: str
    verified: bool
    worker_id: UUID
    size: int
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    ordered_interface_hash: int
    name: str
    description: str | None
    meta_data: str | None
    created: datetime
    updated: datetime
    preload_function: Callable[[], None]
    fitness_function: FitnessFunction
    survivability_function: SurvivabilityFunction
    preload_import: str
    fitness_import: str
    survivability_import: str


class PopulationsConfig(TypedDict):
    """Type definition."""
    configs: NotRequired[list[PopulationConfig]]
    error_on_commit_hash_mismatch: NotRequired[bool]


class PopulationsConfigNorm(TypedDict):
    """Type definition."""
    configs: list[PopulationConfigNorm]
    error_on_commit_hash_mismatch: bool
