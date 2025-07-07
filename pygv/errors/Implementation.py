"""Implementation error classes for PyGV."""


class UnimplementedTransformation(Exception):
    """Raised when a transformation is not implemented."""

    def __init__(self, trans, message="Transformation is not supported"):
        self.message = f"{message} ({trans})"
        super().__init__(self.message)


class UnimplementedBinStat(Exception):
    """Raised when a bin statistic is not implemented."""

    pass
