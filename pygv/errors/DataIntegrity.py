"""Data integrity error classes for PyGV."""


class ChromosomeNotInReference(Exception):
    """Raised when a chromosome is not found in the reference."""

    pass


class BamIndexDoesntExists(Exception):
    """Raised when a BAM index file doesn't exist."""

    pass


class InvaildRegion(Exception):
    """Raised when an invalid region is specified."""

    pass
