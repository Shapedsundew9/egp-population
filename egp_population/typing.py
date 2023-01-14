"""Population Types"""
from typing import TypedDict, Callable, Iterable, Any, NotRequired
from uuid import UUID
from datetime import datetime
from egp_types.xGC import xGC


class Population(TypedDict):
    uid: NotRequired[int]
    worker_id: NotRequired[UUID]
    size: int
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    characterize: Callable[[xGC], tuple[float, float]]
    recharacterize: Callable[[Iterable[xGC]], list[float] | tuple[float, ...]]
    oih: NotRequired[int]
    name: str
    description: NotRequired[str]
    meta_data: NotRequired[str]
    created: NotRequired[datetime]


class PopulationNorm(TypedDict):
    uid: int
    worker_id: UUID
    fitness_function_hash: bytes
    size: int
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    characterize: NotRequired[Callable[[xGC], Any]]
    recharacterize: NotRequired[Callable[[Iterable[xGC]], list[Any] | tuple[Any]]]
    oih: int
    name: str
    description: str
    meta_data: str
    created: datetime


Populations = dict[str, Population]
PopulationsNorm = dict[str, PopulationNorm]
