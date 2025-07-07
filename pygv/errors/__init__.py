"""Error classes for PyGV."""

from .DataIntegrity import (
    BamIndexDoesntExists,
    ChromosomeNotInReference,
    InvaildRegion,
)
from .Formatter import NonStandardBed
from .Implementation import (
    UnimplementedBinStat,
    UnimplementedTransformation,
)

__all__ = [
    "ChromosomeNotInReference",
    "BamIndexDoesntExists",
    "InvaildRegion",
    "NonStandardBed",
    "UnimplementedTransformation",
    "UnimplementedBinStat",
]
