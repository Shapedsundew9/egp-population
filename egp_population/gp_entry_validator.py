"""Validate & normalise JSON Genetic Code definitions."""

from copy import deepcopy
from json import load
from os.path import dirname, join
from cerberus import TypeDefinition
from egp_genomic_library.gms_entry_validator import _gms_entry_validator, GMS_ENTRY_SCHEMA, merge
from egp_physics.gc_graph import gc_graph


GP_ENTRY_SCHEMA = deepcopy(GMS_ENTRY_SCHEMA)
with open(join(dirname(__file__), "formats/gp_entry_format.json"), "r") as file_ptr:
    merge(GP_ENTRY_SCHEMA, load(file_ptr))


class _gp_entry_validator(_gms_entry_validator):

    types_mapping = _gms_entry_validator.types_mapping.copy()
    types_mapping['igraph'] = TypeDefinition('igraph', (gc_graph,), ())

    # TODO: Make errors ValidationError types for full disclosure
    # https://docs.python-cerberus.org/en/stable/customize.html#validator-error

    def _check_with_valid_ancestor_a_ref(self, field, value):
        if value is None and self.document['generation']:
            self._error(field, f'GC has no primary parent (ancestor A) but is not a codon (0th generation).')
        if value is not None and not self.document['generation']:
            self._error(field, f'GC has a primary parent (ancestor A) but is a codon (0th generation).')

    def _check_with_valid_ancestor_b_ref(self, field, value):
        if value is not None and self.document['ancestor_a_ref'] is None:
            self._error(field, f'GC has a secondary parent (ancestor B) but no primary parent (ancestor A).')

    def _check_with_valid_gca_ref(self, field, value):
        if 'A' in self.document['graph'] and value is None:
            self._error(field, f'graph references row A but gca_ref is None.')

    def _check_with_valid_gcb_ref(self, field, value):
        if 'B' in self.document['graph'] and value is None:
            self._error(field, f'graph references row B but gcb_ref is None.')
        if value is not None and self.document['gca'] is None:
            self._error(field, f'gcb_ref is defined but gca_ref is None.')

    def _check_with_valid_pgc_ref(self, field, value):
        if self.document['generation'] and value is None:
            self._error(field, f'Generation is > 0 but pgc_ref is None.')
        if not self.document['generation'] and value is not None:
            self._error(field, f'Generation is 0 but pgc_ref is defined as {value}.')
        if self.document['ancestor_a'] is None and value is not None:
            self._error(field, f'GC has no primary parent (ancestor A) but pgc_ref is defined as {value}.')


gp_entry_validator = _gp_entry_validator(GP_ENTRY_SCHEMA)
