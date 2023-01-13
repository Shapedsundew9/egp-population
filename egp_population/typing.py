"""Population Types"""
from typing import TypedDict
from uuid import UUID
from datetime import datetime
from egp_types.ep_type import vtype


class Population(TypedDict):
    uid: int
    worker_id: UUID
    size: int
    inputs: list[str]
    outputs: list[str]
    oih: bytes
    name: str
    description: str
    meta_data: str
    created: datetime
    vt: vtype
    

Populations = dict[int, Population]