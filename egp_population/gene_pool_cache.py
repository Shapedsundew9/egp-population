"""Gene Pool Cache.

The gene pool cache is a space and time optimised store of GC's. It is designed to
be multi-process friendly.

Naively, the gene pool cache could be implemented as a dictionary with reference keys.
This would be fast but does not scale well. Python dictionaries use huge amounts
of memory and are updated in a spatially broad manner requiring subprocesses to maintain
an almost full copy even if most entries are only read.

The gene pool cache as implemented here maintains a dictionary like interface but takes
advantage of some GC structural design choices to efficiently store data in numpy arrays and
fast index them with a trivial (at least as fast as a dictionary) lookup. It also takes
advantage of GC usgae patterns to cluster stable and volatile GC data which makes efficient
use of OS CoW behaviour.
"""

from numpy import zeros, int16, int32, int64, bool_
from logging import DEBUG, NullHandler, getLogger


# Logging
_logger = getLogger(__name__)
_logger.addHandler(NullHandler())
_LOG_DEBUG = _logger.isEnabledFor(DEBUG)


# Critical constants
# SP_UID_MSB: Sub-Process Unique IDentifier Most Significant Bit
# GPC_UID_MSB: Gene Pool Class Unique IDentifier Most Significant Bit
# GCI_UID_MSB: Genetic Code Index Unique IDentifier Most Significant Bit
_SP_UID_MSB = 63
_SP_UID_LSB = 32
_GPC_UID_MSB = 31
_GPC_UID_LSB = 30
_GCI_UID_MSB = 29
_GCI_UID_LSB = 0
STABLE_CLASS = 0
VOLATILE_CLASS = 1
PGC_CLASS = 2
POPULATION_CLASS = 3

# Masks for top level reference bit fields
_SPB_MASK = ((1 << (_SP_UID_MSB + 1)) - 1) ^ ((1 << _SP_UID_LSB) - 1)
_GPCB_MASK = ((1 << _SP_UID_LSB) - 1) ^ ((1 << _GPC_UID_LSB) - 1)
_GCIB_MASK = (1 << _GPC_UID_LSB) - 1

# After using the _GPCB_MASK on a reference these values can be matched to determine the class
_STABLE_CLASS = STABLE_CLASS << _GPC_UID_LSB
_VOLATILE_CLASS = VOLATILE_CLASS << _GPC_UID_LSB
_PGC_CLASS = PGC_CLASS << _GPC_UID_LSB
_POPULATION_CLASS = POPULATION_CLASS << _GPC_UID_LSB
_CLASSES = (
    _STABLE_CLASS,
    _VOLATILE_CLASS,
    _PGC_CLASS,
    _POPULATION_CLASS
)

# Number of non-population indicies
_NUM_INDICES = 1 << _GPC_UID_LSB

# Population class specific decode of the GP Class index
_POPULATION_UID_WIDTH = 8
_POPULATION_UID_LSB = 22
_POPULATION_NUM = 1 << _POPULATION_UID_WIDTH
_POPULATION_UID_MASK = ((1 << _POPULATION_UID_WIDTH) - 1) << _GPC_UID_LSB
_POPULATION_INDEX_MASK = _GCIB_MASK ^ _POPULATION_UID_MASK

# Used in the Gene Pool to bound the sub-process UIDs
_SPB_UID_WIDTH = _SP_UID_MSB - _GPC_UID_MSB
SP_UID_LOWER_LIMIT = -(1 << (_SPB_UID_WIDTH - 1))
SP_UID_UPPER_LIMIT = (1 << (_SPB_UID_WIDTH - 1)) - 1

# Check the masks
_logger.debug(f'Sub-process UID mask: {_SPB_MASK:016X}')
_logger.debug(f'GP class mask: {_GPCB_MASK:016X}')
_logger.debug(f'GP index mask: {_GCIB_MASK:016X}')
_logger.debug(f'Population UID mask: {_POPULATION_UID_MASK:016X}')
_logger.debug(f'Population index mask: {_POPULATION_INDEX_MASK:016X}')

# FIXME: This is temporary
# The Gene Pool Cache (GPC) has a dictionary like interface plus some additional
# access functions for better performance of common operations. An actual dictionary
# uses way too many resources but is easier to implement.
# In the short term (and may be in the long term as a validation reference)
# the GPC is implemented as a dictionary and the 'optimised' access
# functions emulated.

# NOTE the GPC is a global object. It is instanciated at the bottom of this file.

# The 'temporary' GPC
class gene_pool_cache(dict):
    pass

# The 'optimised' GPC
class _gpc():
    """Mapping to a GC stored in the gene pool cache."""

    __slots__ = '_gp', '_ref'

    def __init__(self, gp, ref):
        # FIXME: Should be able to get rid of the _gp reference.
        self._gp = gp
        self._ref = ref

    def __getitem__(self, k):
        """Return the value for key k.

        Args
        ----
        k (str): A valid GPC field.

        Returns
        -------
        (object): Value of field.
        """
        return self._gp.get_field(self._ref, k)

    def __setitem__(self, k, v):
        """Set the value for key k.

        Args
        ----
        k (str): A valid GPC field.
        v (object): The correct type object to set for field k.
        """
        self._gp.set_field(self._ref, k, v)

    def items(self):
        """Iterate through key:value pairs."""
        return [(f, self._gp.get_field(self._ref, f)) for f in self._gp._fields]

    def keys(self):
        """Iterate through keys."""
        return self._gp._fields

    def values(self):
        """Iterate through values."""
        return [self._gp.get_field(self._ref, f) for f in self._gp._fields]


class _class_data():

    def __init__(self, config, size_increment=4096):
        """Initialise the GP Class data structure.

        Args
        ----
        config(iter(dict)): The pypgtable table config.
        size_increment(int): The number of rows to add to the data storage
                       when reallocation is needed. If a row is a ~1kB then
                       this is ~4MB a pop.
        """
        self._keys = {}
        type_map = []
        type_columns = {}
        self.size_increment = size_increment
        for column in config:
            if column['type'] not in type_map:
                type_map.append(column['type'].upper())
            type_columns[column['type']] = type_columns.get[column['type'], 0] + 1
            self._keys[column['name']] = (type_map.index(column['type']), type_columns[column['type']])
            _logger.debug(f'Field {column["name"]} type {column["type"]} at _type_data'
                          f'[{type_map.index(column["type"])}][{type_columns[column["type"]]}]')
        self._type_data = [self._data_object(typ, type_columns[type]) for typ in type_map]


    def __getitem__(self, field):
        type_idx, field_idx = self._keys[field]
        return self._type_data[type_idx][field_idx]


    def _data_object(self, typ, width):
        """Create a data object for the type typ with width columns.

        Args
        ----
        typ(str): A valid pypgtable column type (postgres data type)
        width(int): The number of column in the data type.

        Returns
        -------
        (object): A data object that meets the criteria for typ
        """
        shape = (self.size_increment,) if width == 1 else (width, self.size_increment)
        if typ == 'SMALLINT' or typ == 'INT2':
            return zeros(shape, dtype=int16)
        if typ == 'INTEGER' or typ == 'INT' or typ == 'INT4':
            return zeros(shape, dtype=int32)
        if typ == 'BIGINT' or typ == 'INT8':
            return zeros(shape, dtype=int64)
        if typ == 'BOOL' or typ == 'BOOLEAN':
            return zeros(shape, dtype=bool_)



# TODO: Implement efficient memory usage
# Current implementation defines the API
class _gene_pool_cache():
    """A high performance local cache for the Gene Pool."""

    def __init__(self, config, size_increment=4096):
        """Initialise the GP Local Cache data structure.

        Args
        ----
        config(iter(dict)): The pypgtable table config.
        size_increment(int): The number of rows to add to the data storage
                       when reallocation is needed. If a row is a ~1kB then
                       this is ~4MB a pop.
        """
        _logger.info('Efficient memory usage NOT YET IMPLEMENTED')

        # Data store in this implementation is a list of dict per class
        self._data = (
            _class_data(config, size_increment),
            _class_data(config, size_increment),
            _class_data(config, size_increment),
            tuple(_class_data(config, size_increment) for _ in range(_POPULATION_NUM))
        )

        # Sub-Process UID for this instance shifted to the right location.
        # i.e. uid_value << _GPCB_UID_MSB. This set by the parent process
        # after forking.
        self.sp_uid = 0

        # The next index for a GC added to a class. When a GC is added
        # to one of these groups and there are no empty indices it must be appended at these positions.
        self._indices = [
            0, # STABLE class
            0, # VOLATILE class
            0, # PGC class
            array((i << _GPC_UID_LSB for i in range(_POPULATION_NUM)), dtype=int64) # POPULATION class
        ]

        # A set of empty indices for every class.
        # Assumption here is the set is small.
        self._empty = (
            None, # STABLE class - nothing can be deleted from here.
            set(), # VOLATILE class
            set(), # PGC class
            tuple(set() for _ in range(_POPULATION_NUM)) # POPULATION class
        )


    def __delitem__(self, k):
        """Mark a GC as deleted.

        The deleted GC's index is added to the empty index set and its

        Args
        ----
        k (int): 64 bit reference key
        """
        gp_class = (k & _GPCB_MASK) >> _GPC_UID_LSB
        empty_set = self._empty[gp_class]
        if gp_class == POPULATION_CLASS:
            empty_set = empty_set[(k & _POPULATION_UID_MASK) >> _POPULATION_UID_LSB]
        if _LOG_DEBUG:
            if gp_class == _STABLE_CLASS:
                raise ValueError(f'Delete a GC from the GP stable class is not permitted.')
            elif gp_class == POPULATION_CLASS:
                population_uid = (k & _POPULATION_UID_MASK) >> _POPULATION_UID_LSB
                _logger.debug(f'GP population {population_uid} empty list length is now {len(empty_set)}.')
            elif gp_class == VOLATILE_CLASS:
                _logger.debug(f'GP volatile sub-GC empty list length is now {len(empty_set)}.')
            else:
                _logger.debug(f'GP PGC empty list length is now {len(empty_set)}.')
        empty_set.add(k & _GCIB_MASK)


    def __getitem__(self, k):
        """Return the GC with reference k.

        Args
        ----
        k (int): 64 bit reference key

        Returns
        -------
        (gcg): GC object
        """
        return gcg(self, k)


    def __setitem__(self, k, v):
        """Set a GC object with reference k from a dict-like GC.

        Args
        ----
        k (int): 64 bit reference key
        v (dict-like): GC as a dict-like object
        """
        if _LOG_DEBUG:
            ref_uid = (k & _SPB_MASK)
            if not ref_uid == self.sp_uid:
                raise ValueError(f'GC reference has an invalid sub-process UID for Gene Pool cache is '
                    f'{(k & _SPB_MASK):016X} and should be {self.sp_uid:016X}')
        # TODO: This can be optimised
        for f, i in v.items():
            self.set_field(k, f, i)


    def __contains__(self, k):
        """Return True if k is in the Gene Pool Cache.

        The look up is sub-process UID agnostic. The assumption is that the
        GP cache is self consistent on creation in the parent process.

        Args
        ----
        k (int): 64 bit reference key

        Returns
        -------
        (bool) True if k in self else False
        """
        gp_class = (k & _GPCB_MASK) >> _GPC_UID_LSB
        if gp_class == POPULATION_CLASS:
            population_uid = (k & _POPULATION_UID_MASK) >> _POPULATION_UID_LSB
            empty_set = self._empty[gp_class][population_uid]
            next_idx = self._indices[gp_class][population_uid]
        else:
            empty_set = self._empty[gp_class]
            next_idx = self._indices[gp_class]

        idx = k & _GCIB_MASK
        return idx < next_idx and idx not in empty_set


    def items(self):
        """Iterate through key:value pairs."""
        for gp_class in range(3):
            max_idx = self._indices[gp_class]
            empty_set = self._empty[gp_class]
            data = self._data[gp_class][0:max_idx]
            for gc in filter(lambda x: x not in empty_set, data):
                yield gc['ref'], gc.values()
        for population_uid in range(_POPULATION_NUM):
            max_idx = self._indices[POPULATION_CLASS][population_uid]
            empty_set = self._empty[POPULATION_CLASS][population_uid]
            data = self._data[POPULATION_CLASS][population_uid]
            for gc in filter(lambda x: x not in empty_set, data):
                yield gc['ref'], gc.values()


    def keys(self):
        """Iterate through keys."""

        # Getting values is expensive so better to implement just for keys.
        for gp_class in range(3):
            max_idx = self._indices[gp_class]
            empty_set = self._empty[gp_class]
            data = self._data[gp_class][0:max_idx]
            for gc in filter(lambda x: x not in empty_set, data):
                yield gc['ref']
        for population_uid in range(_POPULATION_NUM):
            max_idx = self._indices[POPULATION_CLASS][population_uid]
            empty_set = self._empty[POPULATION_CLASS][population_uid]
            data = self._data[POPULATION_CLASS][population_uid]
            for gc in filter(lambda x: x not in empty_set, data):
                yield gc['ref']


    def values(self):
        """Iterate through values."""

        # Values are expensive to get the overhead for keys is small.
        for k, v in self.items():
            yield v


    def set_sp_uid(self, spuid):
        """Define the Sub-Process UID.

        Args
        ----
        uid (int): Sub-Process UId
        """
        _logger.info(f'Sub-process UID set to {spuid:08X}')
        self.sp_uid = (spuid << _SP_UID_LSB) & _SPB_MASK


    def ref(self, gp_class=VOLATILE_CLASS, population_uid=0):
        """Return a new reference.

        Args
        ----
        gp_class (int): one of STABLE_CLASS, VOLATILE_CLASS, PGC_CLASS, POPULATION_CLASS
        population(int): A valid population UID
        sgc (bool): True if the reference is for a sub-GC otherwise it is for a PGC

        Returns
        -------
        (int) A 64 bit reference
        """
        if gp_class == POPULATION_CLASS:
            empty_set = self._empty[gp_class][population_uid]
            indices = self._indices[gp_class]
            pos = population_uid
        else:
            empty_set = self._empty[gp_class]
            indices = self._indices
            pos = gp_class

        if empty_set:
            idx = empty_set.pop()
        else:
            idx = indices[pos]
            indices[pos] += 1

        return self.sp_uid | _CLASSES[gp_class] | idx


    def get(self, k, default=None):
        """Return a GC dict-like object if it exists else return default.

        Args
        ----
        k (int): 64 bit reference key
        default (object): If k not in self then return this object

        Returns
        -------
        (dict-like): GC object
        """
        return self[k] if k in self else default


    def setdefault(self, k, default=None):
        """If key k is not in self, insert key k with value default.

        Args
        ----
        k (int): 64 bit reference key
        default (object): Object to put in dictionary

        Returns
        -------
        (object): value of key k
        """
        if k not in self:
            self[k] = default
        return self[k]


    def update(self, x):
        """Update self with dict-like of dict-like objects x.

        Args
        ----
        x (dict-like(dict-like)): Objects to update self with.
        """
        for k, v in x.items():
            self[k] = v


    def _find_field(self, k, f):
        """Find the data object and index for reference k, field f.

        Args
        ----
        k (int): Reference
        f (str): Field name

        Returns
        -------
        (object, int): An indexable data object and the index into it where
                       the object k, f references is.
        """
        gp_class = (k & _GPCB_MASK) >> _GPC_UID_LSB
        if gp_class == POPULATION_CLASS:
            population_uid = (k & _POPULATION_UID_MASK) >> _POPULATION_UID_LSB
            data = self._data[POPULATION_CLASS][population_uid]
            idx = k & _POPULATION_INDEX_MASK
        else:
            data = self._data[gp_class]
            idx = k & _GCIB_MASK
        return data[f], idx


    def get_field(self, k, f):
        """Return the value of self[k][f].

        Args
        ----
        k (int): 64 bit reference key
        f (int): Field index

        Returns
        -------
        (object): self[k][f]
        """
        data, idx = self._find_field(k, f)
        return data[idx]


    def set_field(self, k, f, v):
        """Set self[k][f] to v.

        Args
        ----
        k (int): 64 bit reference key
        f (int): Field index
        default (object): Object to put in dictionary

        Returns
        -------
        (object): value of k[f]
        """
        data, idx = self._find_field(k, f)
        data[idx] = v


    def getdefault_field(self, k, f, default=None):
        """Return the value of self[k][f] if it exists else return default.

        default is returned if either k or k[f] do not exist.

        Args
        ----
        k (int): 64 bit reference key
        f (int): Field index
        default (object): If self[k][f] not in self then return this object

        Returns
        -------
        (object): self[k][f] or default
        """
        if k in self:
            return self.get_field(k, f)
        return default


    def get_field_array(self, gp_class, field, population_uid=0):
        """Return an iterable of the field f for GP GC class gp_class.

        If gp_class == POPULATION_CLASS then the array of f is for population p.
        The order of values returned is guaranteed to be the same between
        calls if self is not modified.

        Args
        ----
        gp_class (int): one of STABLE_CLASS, VOLATILE_CLASS, PGC_CLASS, POPULATION_CLASS
        field (str): Field name
        population_uid (int): Population UID

        Returns
        -------
        iter(object): Values of field for all GCs in the class / population
        """
        if _LOG_DEBUG:
            if gp_class == POPULATION_CLASS and self._empty[gp_class][population_uid]:
                raise ValueError(f'Population UID {population_uid} has empty values. Field array {field} is not contiguous!')
            elif self._empty[gp_class]:
                raise ValueError(f'Class {gp_class} has empty values. Field array {field} is not contiguous!')

        if gp_class == POPULATION_CLASS:
            return self._data[POPULATION_CLASS][population_uid][field][0:self._indices[POPULATION_CLASS][population_uid]]
        return self._data[gp_class][field][0:self._indices[gp_class]]
