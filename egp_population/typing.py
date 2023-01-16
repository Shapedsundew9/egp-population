"""Population Types"""
from typing import TypedDict, Callable, Iterable, NotRequired
from uuid import UUID
from datetime import datetime
from egp_types.xGC import xGC
from egp_types.ep_type import vtype


class Population(TypedDict):
    uid: NotRequired[int]
    worker_id: NotRequired[UUID]
    size: int
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    fitness_function: Callable[[xGC], float]
    survivability_function: Callable[[Iterable[xGC]], list[float] | tuple[float, ...]]
    ordered_interface_hash: NotRequired[int]
    name: str
    description: NotRequired[str | None] 
    meta_data: NotRequired[str | None]
    vt: NotRequired[vtype]


class PopulationNorm(TypedDict):
    uid: int
    worker_id: UUID
    fitness_function_hash: bytes
    size: int
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    fitness_function: Callable[[xGC], float]
    survivability_function: Callable[[Iterable[xGC]], list[float] | tuple[float]]
    ordered_interface_hash: int
    name: str
    description: str | None
    meta_data: str | None
    created: datetime
    vt: vtype


Populations = dict[str, Population]
PopulationsNorm = dict[str, PopulationNorm]
