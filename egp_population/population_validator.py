"""Validate & normalise JSON Genetic Code definitions."""

from json import load
from os.path import dirname, join
from typing import Any

from egp_utils.base_validator import base_validator
from egp_types.ep_type import validate, ordered_interface_hash, interface_definition, vtype


with open(join(dirname(__file__), "formats/gp_entry_format.json"), "r", encoding="utf8") as file_ptr:
    POPULATION_ENTRY_SCHEMA: dict[str, dict[str, Any]] = load(file_ptr)


class _population_entry_validator(base_validator):

    def _check_with_valid_type(self, field: str, value: Any) -> None:
        vt: vtype = self.document.get('vt', vtype.EP_TYPE_STR)
        if not validate(value, vt):
            self._error(field, f'ep_type {value} does not exist with vtype {vt}.')

    def _normalize_default_setter_set_oih(self, document) -> int:
        vt: vtype = document.get('vt', vtype.EP_TYPE_STR)
        o_def: tuple[tuple[int, ...], list[int], bytes] = interface_definition(self.document['outputs'], vt)
        i_def: tuple[tuple[int, ...], list[int], bytes] = interface_definition(self.document['inputs'], vt)
        return ordered_interface_hash(i_def[1], o_def[1], i_def[2], o_def[2])


population_entry_validator: _population_entry_validator = _population_entry_validator(POPULATION_ENTRY_SCHEMA)
