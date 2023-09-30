"""Validate & normalise JSON Genetic Code definitions."""

from json import load
from os.path import dirname, join
from typing import Any
from datetime import datetime
from egp_utils.base_validator import base_validator
from egp_types.ep_type import validate, vtype


with open(join(dirname(__file__), "formats/population_entry_format.json"), "r", encoding="utf8") as file_ptr:
    POPULATION_ENTRY_SCHEMA: dict[str, dict[str, Any]] = load(file_ptr)


class _population_entry_validator(base_validator):
    def _check_with_valid_ep_type(self, field: str, value: str) -> None:
        val_t: vtype = self.document.get("vt", vtype.EP_TYPE_STR)
        if not validate(value, val_t):
            self._error(field, f"ep_type {value} does not exist with vtype {val_t}.")

    def _check_with_valid_created(self, field: str, value: datetime) -> None:
        if value > datetime.utcnow():
            self._error(
                field,
                "Created date-time cannot be in the future. Is the system clock correct?",
            )
        if self.document.get("updated") is not None:
            if self.document["updated"] < value:
                self._error(field, "A record cannot be updated before it has been created.")

    def _check_with_valid_updated(self, field: str, value: datetime) -> None:
        if value > datetime.utcnow():
            self._error(
                field,
                "Updated date-time cannot be in the future. Is the system clock correct?",
            )
        if value < self.document["created"]:
            self._error(field, "A record cannot be updated before it has been created.")


population_entry_validator: _population_entry_validator = _population_entry_validator(POPULATION_ENTRY_SCHEMA)
