"""Population Types"""
from typing import TypedDict, NotRequired, Iterable, Callable
from uuid import UUID
from datetime import datetime
from egp_types.xGC import xGC
from egp_types.egp_typing import FitnessFunction, SurvivabilityFunction


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


# Used in the 'cast' of PopulationConfig to PopulationConfigNorm
# Seems like a lot of hoops to jump through to get type checking.
def _mock_plf() -> None:
    pass


def _mock_ff(dummy: Iterable[xGC]) -> None:
    """Function is never executed."""
    for dumb in dummy:
        dumb[''] = None


def _mock_sf(dummy1: Iterable[xGC], dummy2: Iterable[xGC]) -> None:
    """Function is never executed."""
    for dumb in dummy1:
        dumb[''] = None
    for dumb in dummy2:
        dumb[''] = None


def cast_pc_to_pcn(population_config: PopulationConfig) -> PopulationConfigNorm:
    """Force a normalized PopulationConfig type to a PopulationConfigNorm."""
    return {
        'uid': population_config.get('uid', int()),
        'population_hash': population_config.get('population_hash', bytes()),
        'git_repo': population_config.get('git_repo', str()),
        'git_url': population_config.get('git_url', str()),
        'git_hash': population_config.get('git_hash', str()),
        'verified': population_config.get('verified', bool()),
        'worker_id': population_config.get('worker_id', UUID()),
        'size': population_config.get('size', int()),
        'inputs': population_config.get('inputs', tuple()),
        'outputs': population_config.get('outputs', tuple()),
        'ordered_interface_hash': population_config.get('ordered_interface_hash', int()),
        'name': population_config.get('name', str()),
        'description': population_config.get('description'),
        'meta_data': population_config.get('meta_data'),
        'created': population_config.get('created', datetime.now()),
        'updated': population_config.get('updated', datetime.now()),
        'preload_function': population_config.get('preload_function', _mock_plf),
        'fitness_function': population_config.get('fitness_function', _mock_ff),
        'survivability_function': population_config.get('survivability_function', _mock_sf),
        'preload_import': population_config.get('preload_import', str()),
        'fitness_import': population_config.get('fitness_import', str()),
        'survivability_import': population_config.get('survivability_import', str())
    }
