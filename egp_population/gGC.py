"""Specialized GC object for the Gene Pool."""

from logging import DEBUG, NullHandler, getLogger
from copy import copy
from egp_physics.gc_type import _GC, interface_definition, is_pgc, _GL_EXCLUDE_COLUMNS, NUM_PGC_LAYERS
from .population_validator import gp_entry_validator
from egp_physics.gc_graph import gc_graph
from egp_physics.execution import remove_callable
from egp_physics.ep_type import vtype


_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)

GPC = None

def set_GPC(gpc):
    global GPC
    GPC = gpc


def gGC(gcs=tuple(), modified=True, population_uid=None, sv=True):
    """Wrap creation of gGC objects.

    This exists to do global scope validation and
    optimisations when a gGC is created. It is a single point of entry
    to the Gene Pool Cache.

    gcs must be self consistent. Any referenced GC must either already be
    in the GPC or in gcs.

    Validations
        1. Does the GC already exist.
        2. Does an allotrope exist. NOT IMPLEMENTED.
        3. Do GCA &| GCB exist in the cache. NOT IMPLEMENTED.

    Optimisations
        1. Consolidate executable code to reduce the size of the call graph. IN PROGRESS.
        2. Identify duplicates. NOT IMPLEMENTED. TODO: Dupes should reduce probability of creation of another dupe.
        3. Alias allotropes. NOT IMPLEMENTED.
    """
    ggcs = []
    for ggc in gcs:
        if _LOG_DEBUG and 'signature' in ggc:
            _logger.debug(f"New gGC with signature: {ggc['signature']} (ref: {_GC.ref_from_sig(ggc['signature'])}) found.")

        if 'ref' in ggc and ggc['ref'] in GPC:
            if _LOG_DEBUG:
                _logger.debug(f"gGC is already in the Gene Pool.")
            ggcs.append(GPC[ggc['ref']])
        else:
            ggc = _gGC(ggc, modified, population_uid, sv)
            GPC[ggc['ref']] = ggc
            ggcs.append(ggc)

    # Sanity
    if _LOG_DEBUG:
        refs = {ggc['ref'] for ggc in ggcs}
        for ggc in ggcs:
            for ref_check in ('gca_ref', 'gcb_ref', 'pgc_ref', 'ancestor_a_ref', 'ancestor_b_ref'):
                if ggc.get(ref_check, None) is not None:
                    if not (ggc[ref_check] in refs or ggc[ref_check] in GPC):
                        assert False, f"GC {ggc['ref']}'s {ref_check} ({ggc[ref_check]}) does not exist!"

    return ggcs

set_next_reference = _GC.set_next_reference


class _gGC(_GC):
    """Gene Pool GC type.

    Gene pool GC types hold a lot of transient data.
    gGC's are stored in the Gene Pool Cache.
    Only 1 instance can exist: Copying or recreating will raise an exception.
    """

    higher_layer_cols = tuple((col for col in filter(lambda x: x[0] == '_', gp_entry_validator.schema.keys())))

    def __init__(self, gc={}, modified=True, population_uid=None, sv=True):
        """Construct.

        Ensure all fields are defined as required for the Gene Pool.
        If gc has a 'signature' then it may have been pulled from the
        genomic library and already exist in the gene pool cache. If so
        the gGC is not recreated.

        Args
        ----
        gc (a _GC dervived object): GC to ensure is eGC compliant.
        sv (bool): Suppress validation. If True the eGC will not be validated on construction.
        """
        super().__init__(gc)
        if _LOG_DEBUG:
            assert self['ref'] is not None
            _logger.debug(f"Added gGC with ref: {self['ref']} to GP cache.")
            _logger.debug(f"GP cache now has {len(GPC)} entries.")
        self['modified'] = modified
        self.setdefault('ac_count', 1)
        self.setdefault('cb')

        # FIXME: Not needed if ref is a function of population_uid.
        self.setdefault('population_uid', population_uid)
        _logger.debug(f"self['population_uid']={self['population_uid']}, population_uid={population_uid}")

        self.setdefault('pgc_ref', self.field_reference('pgc'))
        self.setdefault('gca_ref', self.field_reference('gca'))
        self.setdefault('gcb_ref', self.field_reference('gcb'))
        self.setdefault('igraph', gc_graph(self['graph']))
        self.setdefault('generation', 0)
        self.setdefault('sms_ref', None)
        self.setdefault('effective_pgc_refs', [])
        self.setdefault('effective_pgc_fitness', [])
        self.setdefault('offspring_count', 0)
        if 'inputs' not in self:
            inputs = self['igraph'].input_if()
            outputs = self['igraph'].output_if()
            _, self['input_types'], self['inputs'] = interface_definition(inputs, vtype.EP_TYPE_INT)
            _, self['output_types'], self['outputs'] = interface_definition(outputs, vtype.EP_TYPE_INT)

        self['num_inputs'] = len(self['inputs'])
        self['num_outputs'] = len(self['outputs'])
        self['callable'] = None

        for col in filter(lambda x: x[1:] in gc.keys(), _gGC.higher_layer_cols):
            gc[col] = copy(gc[col[1:]])

        # PGCs have special fields in the Gene Pool
        if is_pgc(self) and 'pgc_f_valid' not in self:
            self['pgc_f_valid'] = [f > 0.0 for f in self['pgc_fitness']]
        else:
            self.setdefault('fitness', 0.0)
            self.setdefault('survivability', 0.0)

        # Remove Genomic Library fields
        for gl_column in filter(lambda x: x in self, _GL_EXCLUDE_COLUMNS):
            del self[gl_column]

        if _LOG_DEBUG:
            # Inefficient to recreate gGC's.
            assert not isinstance(gc, _gGC)

            # Avoid a copy of gGC which causes issues with __del__() of the 'exec' function.
            gc = dict(self)
            _logger.debug(f'Validating gGC: {self} using dictionary: {gc}')
            _logger.debug(','.join([f'{k}: {type(v)}({v})' for k, v in gc.items()]))
            assert gp_entry_validator(gc), f'gGC invalid:\n{gp_entry_validator.error_str()}.'
            assert self['ref'] != self['gca_ref'], 'GC ref == GCA ref. A GC cannot self reference.'
            assert self['ref'] != self['gcb_ref'], 'GC ref == GCB ref. A GC cannot self reference.'

    def __del__(self):
        """Make sure the 'exec' function is cleaned up."""
        if _LOG_DEBUG:
            _logger.debug(f"Deleting {self['ref']}.")
        remove_callable(self)

    def __copy__(self):
        """Make sure we do not copy gGCs."""
        assert False, f"Shallow copy of gGC {self['ref']}"

    def __deepcopy__(self):
        """Make sure we do not copy gGCs."""
        assert False, f"Deep copy of gGC {self['ref']}"
