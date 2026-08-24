"""Shared types and validators for track parameters."""

from __future__ import annotations

from typing import Annotated, Any, Callable, Literal, Optional

from matplotlib.colors import is_color_like
from pydantic import AfterValidator, Field

ShowMode = Literal["expanded", "collapsed"]
PlotType = Literal["line", "bar"]
StackOrder = Literal["big_on_top", "small_on_top", "fixed"]
ColorReadsBy = Literal[
    "insert size",
    "pair orientation",
    "insert size and pair orientation",
    "read strand",
    "first of pair strand",
    "read group",
    "sample",
    "library",
    "movie",
    "ZMW",
    "tag",
    "no color",
]

BIN_STATS = frozenset(
    {"mean", "std", "median", "count", "sum", "min", "max"}
)


def _validate_color(value: Any) -> Any:
    if value is None:
        return value
    if not is_color_like(value):
        raise ValueError(f"Invalid color value: {value}")
    return value


def _validate_color_sequence(value: Any) -> Any:
    if value is None:
        return value
    if isinstance(value, (str, bytes)) or not hasattr(value, "__iter__"):
        raise ValueError("Expected a sequence of colors")
    validated = tuple(_validate_color(item) for item in value)
    if len(validated) == 0:
        raise ValueError("Expected a non-empty sequence of colors")
    return validated


def _validate_filter(value: Any) -> Any:
    if value is None or callable(value):
        return value
    raise ValueError("filters must be None or a callable")


Color = Annotated[Any, AfterValidator(_validate_color)]
ColorSequence = Annotated[Any, AfterValidator(_validate_color_sequence)]
Alpha = Annotated[float, Field(ge=0, le=1)]
PositiveFloat = Annotated[float, Field(gt=0)]
FilterFn = Annotated[Optional[Callable[..., Any]], AfterValidator(_validate_filter)]
