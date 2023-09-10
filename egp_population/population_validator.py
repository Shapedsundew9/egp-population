"""Validate & normalise JSON Genetic Code definitions."""

from json import load
from os.path import dirname, join
from typing import Any, LiteralString
from hashlib import sha256
from datetime import datetime
from egp_utils.base_validator import base_validator
from egp_types.ep_type import (
    validate,
    ordered_interface_hash,
    interface_definition,
    vtype,
)


with open(
    join(dirname(__file__), "formats/gp_entry_format.json"), "r", encoding="utf8"
) as file_ptr:
    POPULATION_ENTRY_SCHEMA: dict[str, dict[str, Any]] = load(file_ptr)


_POPULATION_HASH_FIELDS: tuple[LiteralString, ...] = (
    "git_repo_url",
    "git_hash",
    "size",
    "ordered_interface_hash",
    "name",
)


class _population_entry_validator(base_validator):
    def _check_with_valid_type(self, field: str, value: str) -> None:
        vt: vtype = self.document.get("vt", vtype.EP_TYPE_STR)
        if not validate(value, vt):
            self._error(field, f"ep_type {value} does not exist with vtype {vt}.")

    def _check_with_valid_git_repo_url(self, field: str, value: str | None) -> None:
        no_url: bool = value is None
        no_hash: bool = self.document.get("git_hash", None) is None
        no_repo: bool = self.document.get("git_repo", None) is None
        if not (no_url == no_hash and no_hash == no_repo):
            self._error(
                field,
                "All of 'git_repo_url', 'git_hash' & 'git_repo' must be defined or all must be None.",
            )

    def _check_with_valid_git_repo(self, field: str, value: bool | None) -> None:
        no_url: bool = self.document.get("git_repo", None) is None
        no_hash: bool = self.document.get("git_hash", None) is None
        no_repo: bool = value is None
        if not (no_url == no_hash and no_hash == no_repo):
            self._error(
                field,
                "All of 'git_repo_url', 'git_hash' & 'git_repo' must be defined or all must be None.",
            )

    def _check_with_valid_git_hash(self, field: str, value: bytes | None) -> None:
        no_url: bool = self.document.get("git_repo_url", None) is None
        no_hash: bool = value is None
        no_repo: bool = self.document.get("git_repo", None) is None
        if not (no_url == no_hash and no_hash == no_repo):
            self._error(
                field,
                "All of 'git_repo_url', 'git_hash' & 'git_repo' must be defined or all must be None.",
            )

    def _check_with_valid_verified(self, field: str, value: bool | None) -> None:
        no_url: bool = self.document.get("git_repo_url", None) is None
        no_hash: bool = self.document.get("git_hash", None) is None
        no_repo: bool = self.document.get("git_repo", None) is None
        if no_url or no_hash or no_repo and value is not None:
            self._error(field, "verified cannot be defined if git_* are not defined.")

    def _check_with_valid_created(self, field: str, value: datetime) -> None:
        if value > datetime.utcnow():
            self._error(
                field,
                "Created date-time cannot be in the future. Is the system clock correct?",
            )
        if self.document.get("updated") is not None:
            if self.document["updated"] < value:
                self._error(
                    field, "A record cannot be updated before it has been created."
                )

    def _check_with_valid_updated(self, field: str, value: datetime) -> None:
        if value > datetime.utcnow():
            self._error(
                field,
                "Updated date-time cannot be in the future. Is the system clock correct?",
            )
        if self.document.get("created") is not None:
            if self.document["updated"] > value:
                self._error(
                    field, "A record cannot be updated before it has been created."
                )

    def _normalize_default_setter_set_oih(self, document) -> int:
        vt: vtype = document.get("vt", vtype.EP_TYPE_STR)
        o_def: tuple[tuple[int, ...], list[int], bytes] = interface_definition(
            self.document["outputs"], vt
        )
        i_def: tuple[tuple[int, ...], list[int], bytes] = interface_definition(
            self.document["inputs"], vt
        )
        return ordered_interface_hash(i_def[1], o_def[1], i_def[2], o_def[2])

    def _normalize_default_setter_set_population_hash(self, document) -> bytes:
        string: str = "".join(
            (str(document.get(field, None)) for field in _POPULATION_HASH_FIELDS)
        )
        return sha256(string.encode()).digest()

    def _normalize_default_setter_set_verified(self, document) -> bool | None:
        no_repo: bool = self.document.get("git_repo_url", None) is None
        no_hash: bool = self.document.get("git_hash", None) is None
        no_256: bool = self.document.get("git_repo", None) is None
        if no_repo or no_hash or no_256:
            return None
        if not no_repo or not no_hash or not no_256:
            return False


population_entry_validator: _population_entry_validator = _population_entry_validator(
    POPULATION_ENTRY_SCHEMA
)
