"""Population Types"""
from typing import TypedDict, NotRequired
from uuid import UUID
from datetime import datetime
from egp_types.typing import FitnessFunction, SurvivabilityFunction


class PopulationConfig(TypedDict):
    uid: NotRequired[int]
    population_hash: NotRequired[bytes]
    git_repo: NotRequired[str]
    git_url: NotRequired[str]
    git_hash: NotRequired[bytes]
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
    fitness_function: NotRequired[FitnessFunction]
    survivability_function: NotRequired[SurvivabilityFunction]
    fitness_import: NotRequired[str]
    survivability_import: NotRequired[str]


class PopulationConfigNorm(TypedDict):
    uid: int
    population_hash: bytes
    git_repo: str
    git_url: str
    git_hash: bytes
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
    fitness_function: FitnessFunction
    survivability_function: SurvivabilityFunction
    fitness_import: str
    survivability_import: str


PopulationsConfig = list[PopulationConfig]
PopulationsConfigNorm = list[PopulationConfigNorm]
