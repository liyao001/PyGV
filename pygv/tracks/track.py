from __future__ import annotations

from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    field_validator,
)

from pygv.errors.DataIntegrity import InvaildRegion
from pygv.tracks.types import Alpha, Color, PositiveFloat

TRACK_MODEL_CONFIG = ConfigDict(
    arbitrary_types_allowed=True,
    validate_assignment=True,
    extra="forbid",
    populate_by_name=True,
)


def _echo(data):
    return data


_TRANSFORMATIONS = {
    "ln": np.log,
    "asinh": np.arcsinh,
    "log2": np.log2,
    "log10": np.log10,
    "log1p": np.log1p,
    "rln": lambda x: -1 * np.log(-1 * x),
    "rlog2": lambda x: -1 * np.log2(-1 * x),
    "rlog10": lambda x: -1 * np.log10(-1 * x),
    "rlog1p": lambda x: -1 * np.log1p(-1 * x),
}

_BIN_STATS = {"mean", "std", "median", "count", "sum", "min", "max"}


class Track(BaseModel):
    """Generic track. Parameters are typed class fields for IDE autocomplete and validation."""

    model_config = TRACK_MODEL_CONFIG

    name: str = Field(default="", kw_only=True, description="Name of the track")
    line_width: float = Field(
        default=1, kw_only=True, description="The default width for lines"
    )
    height: PositiveFloat = Field(
        default=1,
        kw_only=True,
        description="Height of the track (unit, relative measurement)",
    )
    color: Color = Field(default="#A1A1A1", kw_only=True, description="Default color")
    edge_color: Color = Field(
        default="#6E6E6E", kw_only=True, description="Edge color"
    )
    alpha: Alpha = Field(default=0.8, kw_only=True, description="Alpha of patches")
    font_color: Color = Field(default="black", kw_only=True, description="Font color")
    font_size: Optional[float] = Field(
        default_factory=lambda: plt.rcParams["font.size"],
        kw_only=True,
        description="Font size",
    )
    y_tick_format: Optional[str] = Field(
        default=None,
        kw_only=True,
        description="String format for ticks on y-axis. For example: `{:.1f}`",
    )
    y_label_rotation: Union[str, float] = Field(
        default="horizontal",
        kw_only=True,
        description="Rotation of y-axis' label",
    )
    y_label_ha: str = Field(
        default="right",
        kw_only=True,
        description="Horizontal alignment about label for y-axis",
    )
    y_label_va: str = Field(
        default="center",
        kw_only=True,
        description="Vertical alignment about label for y-axis",
    )
    hide_left_spine: bool = Field(
        default=False,
        kw_only=True,
        description="Hide the left spine of the track axis",
    )
    inward_yticks: bool = Field(
        default=False,
        kw_only=True,
        validation_alias=AliasChoices("inward_yticks", "inward_ticks"),
        description=(
            "Plot y-ticks strictly inside each track. "
            "To apply this to all tracks, set inward_ticks=True on GenomeViewer."
        ),
    )

    _ax: Any = PrivateAttr(default=None)
    _highlight_starts: list = PrivateAttr(default_factory=list)
    _highlight_ends: list = PrivateAttr(default_factory=list)
    _highlight_colors: list = PrivateAttr(default_factory=list)
    _highlight_alphas: list = PrivateAttr(default_factory=list)

    def layout_height(self) -> float:
        """Height used by GenomeViewer when allocating subplot space."""
        return self.height

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        """Called before :func:`~pygv.tracks.track.Track._draw_track`."""
        if self.font_size is None:
            self.font_size = plt.rcParams["font.size"]
        inward = kwargs.get("inward_ticks", kwargs.get("inward_yticks", None))
        if inward is not None:
            self.inward_yticks = inward

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """Draw track chrome (axis limits, spines, coordinate labels)."""
        self._ax = ax
        n_ticks = kwargs.get("n_ticks", None)
        hide_coords = kwargs.get("hide_coordinates", False)
        hide_chr_name = kwargs.get("hide_chromosome_name", None)
        if start <= end:
            self._ax.set_xlim((start, end))
        else:
            raise InvaildRegion("start of the region must be smaller than the end")
        if index != 0:
            self._ax.spines["top"].set_visible(False)
            self._ax.set_xticks([])
            self._ax.set_xticks([], minor=True)
        elif hide_coords:
            self._ax.spines["top"].set_visible(False)
            self._ax.tick_params(
                axis="x",
                which="both",
                top=False,
                bottom=False,
                labeltop=False,
                labelbottom=False,
            )
        else:
            self._ax.xaxis.set_ticks_position("top")
            self._ax.spines["top"].set_position(("outward", 10))
            self._ax.spines["top"].set_linewidth(2)
            if n_ticks is not None:
                ticks = np.linspace(start, end, n_ticks, dtype=int)
            else:
                ticks = [t for t in self._ax.get_xticks() if start <= t <= end]

            if ticks[-1] - ticks[1] <= 1e3:
                labels = [f"{x:,.0f}" for x in ticks]
                labels[-1] += " bp"
            elif ticks[-1] - ticks[1] <= 4e5:
                labels = [f"{x / 1000.0:,.0f}" for x in ticks]
                labels[-1] += " Kb"
            else:
                labels = [f"{x / 1000000.0:,.1f} " for x in ticks]
                labels[-1] += " Mbp"

            if not hide_chr_name:
                self._ax.set_title(chromosome)
            self._ax.set_xticks(ticks)
            self._ax.set_xticklabels(labels)

        self._ax.spines["bottom"].set_visible(False)
        self._ax.spines["right"].set_visible(False)

        if self.hide_left_spine:
            self._ax.spines["left"].set_visible(False)

        if self.name is not None:
            self._ax.set_ylabel(
                self.name,
                rotation=self.y_label_rotation,
                ha=self.y_label_ha,
                va=self.y_label_va,
            )

    def set_highlight_regions(self, starts, ends, colors=(), alpha_vals=()):
        """Set highlight regions."""
        if isinstance(starts, int) or isinstance(starts, float):
            starts = [starts]
        if isinstance(ends, int) or isinstance(ends, float):
            ends = [ends]
        n_starts = len(starts)
        n_ends = len(ends)
        if n_starts != n_ends:
            raise ValueError(
                "The number of start positions must be equal to the number of end positions."
            )
        if isinstance(colors, str):
            colors = (colors,)
        if not all(map(is_color_like, colors)):
            raise ValueError("invalid color value(s).")
        n_colors = len(colors)
        if n_colors == 0:
            self._highlight_colors = ("#FFFC66",) * n_starts
        elif n_colors == 1:
            self._highlight_colors = colors * n_starts
        elif n_colors == n_starts:
            self._highlight_colors = colors
        else:
            raise ValueError(
                "You need to provide either no, one, or colors for each region."
            )
        n_alpha = len(alpha_vals)
        if n_alpha == 0:
            self._highlight_alphas = (0.5,) * n_starts
        elif n_alpha == 1:
            self._highlight_alphas = alpha_vals * n_starts
        elif n_alpha == n_starts:
            self._highlight_alphas = alpha_vals
        else:
            raise ValueError(
                "You need to provide either no, one, or alpha values for each region."
            )
        self._highlight_starts = starts
        self._highlight_ends = ends

    def add_highlight_region(self, start, end):
        """Add one highlight region."""
        self._highlight_starts.append(float(start))
        self._highlight_ends.append(float(end))
        colors = list(self._highlight_colors)
        alphas = list(self._highlight_alphas)
        colors.append("#FFFC66")
        alphas.append(0.5)
        self._highlight_colors = colors
        self._highlight_alphas = alphas

    def remove_highlight(self):
        """Remove highlight zones."""
        self._highlight_starts = []
        self._highlight_ends = []
        self._highlight_colors = []
        self._highlight_alphas = []

    def _post_plot_hook(self, chromosome, start, end, ax, index=1, **kwargs):
        """Draw highlight spans after the track is drawn."""
        n_starts = len(self._highlight_starts)
        n_ends = len(self._highlight_ends)
        if n_starts == n_ends and n_starts > 0:
            for s, e, c, a in zip(
                self._highlight_starts,
                self._highlight_ends,
                self._highlight_colors,
                self._highlight_alphas,
            ):
                self._ax.axvspan(s, e, color=c, alpha=a, linewidth=0, zorder=-1)


class AnnotationTrack(Track):
    """Annotation track for interval features such as genes and peaks."""

    track: str = Field(description="Path to the annotation file")
    patch_height: float = Field(
        default=1, kw_only=True, description="Height of patches (for exons/blocks)"
    )
    allowed_feature_lanes: Optional[int] = Field(
        default=None,
        kw_only=True,
        description=(
            "Max amount of feature lanes to be plotted. For example, if a region has "
            "12 overlapping features, they will be plotted into 12 separate lanes. "
            "If you set this to 2, only two lanes will be shown."
        ),
    )
    font_box_alpha: Alpha = Field(
        default=0.75,
        kw_only=True,
        description="Transparent/alpha for text boxes labeling gene names",
    )
    lane_space: float = Field(
        default=0.25, kw_only=True, description="Extra spaces between lanes"
    )
    features_per_lane: int = Field(
        default=3, kw_only=True, description="Features per lane"
    )
    line_color: Color = Field(
        default="black", kw_only=True, description="Line color"
    )
    arrow_interval: float = Field(
        default=5, kw_only=True, description="Intervals between arrows"
    )
    padding_left: float = Field(
        default=0,
        kw_only=True,
        description=(
            "Extra left padding so feature names do not overlap. An integer is a "
            "base-pair distance; a float between 0 and 1 is a fraction of the window."
        ),
    )
    padding_right: float = Field(
        default=0,
        kw_only=True,
        description="Extra right padding for feature labels",
    )
    show_name: bool = Field(
        default=True,
        kw_only=True,
        description="Print names of genomic regions if available",
    )
    hide_visual_dup: bool = Field(
        default=False,
        kw_only=True,
        description=(
            "Hide features which are duplicates of other features in the current window"
        ),
    )

    _plot_thickness: int = PrivateAttr(default=0)
    _plot_block: int = PrivateAttr(default=0)
    _small_relative: float = PrivateAttr(default=0)
    _lane_registries: list = PrivateAttr(default_factory=list)

    def __init__(self, track: str, **data: Any) -> None:
        super().__init__(track=track, **data)

    def layout_height(self) -> float:
        return max(len(self._lane_registries), 1) * self.height

    def _plot_gene_direction(self, ax, xpos, ypos, strand, **kwargs):
        """Draw a broken line indicating strand: `>` for `+`, `<` for `-`."""
        if strand == ".":
            return

        if strand == "+":
            xdata = [
                xpos - self._small_relative / 3,
                xpos + self._small_relative / 3,
                xpos - self._small_relative / 3,
            ]
        else:
            xdata = [
                xpos + self._small_relative / 3,
                xpos - self._small_relative / 3,
                xpos + self._small_relative / 3,
            ]

        ydata = [ypos - 1 / 5, ypos, ypos + 1 / 5]
        ax.add_line(
            Line2D(xdata, ydata, color=self.line_color, linewidth=self.line_width)
        )


class NumericalTrack(Track):
    """Numerical track for continuous signals."""

    min_val: Optional[float] = Field(
        default=None,
        kw_only=True,
        description="Minimum value to be plotted. By default, all signals are plotted.",
    )
    max_val: Optional[float] = Field(
        default=None,
        kw_only=True,
        description="Maximum value to be plotted. By default, all signals are plotted.",
    )
    show_range: bool = Field(
        default=True, kw_only=True, description="Whether to show y-axis range ticks"
    )
    n_bins: Optional[int] = Field(
        default=None,
        kw_only=True,
        description=(
            "Number of bins to apply. If a positive number is set, the window will "
            "be separated in bins and stat_method will be applied."
        ),
    )
    stat_method: Optional[str] = Field(
        default=None,
        kw_only=True,
        description="Statistical method for binning windows",
    )
    data_transform: Any = Field(
        default=None,
        kw_only=True,
        validate_default=True,
        validation_alias=AliasChoices("data_transform", "transformation"),
        description=(
            "Function for data transformation. None, a named transform "
            "(`asinh`, `ln`, `log2`, `log10`, `log1p`, or `r`-prefixed variants), "
            "or a callable."
        ),
    )
    convert_nan_to_num: Any = Field(
        default=np.nan_to_num,
        kw_only=True,
        validate_default=True,
        description="Function mapping NaN values. None leaves values unchanged.",
    )
    scale: float = Field(
        default=1,
        kw_only=True,
        description="Normalization factor for signals (e.g. RPM)",
    )
    label_masked_peak: bool = Field(
        default=True,
        kw_only=True,
        description="Whether to label capped overflow signals",
    )
    overflow_label_format: Optional[str] = Field(
        default="{:.1f}",
        kw_only=True,
        description="String format for labeling overflow loci",
    )
    overflow_label_auto_adjust: bool = Field(
        default=False,
        kw_only=True,
        description="Automatically place text labels for overflow signals",
    )
    equal_space_for_pos_neg_ranges: bool = Field(
        default=False,
        kw_only=True,
        validation_alias=AliasChoices(
            "equal_space_for_pos_neg_ranges", "draw_y_independently"
        ),
        description="Force positive and negative y-ranges to occupy equal space",
    )
    skip_label_for_zero: bool = Field(default=False, kw_only=True)

    _yscale_func: Any = PrivateAttr(default=None)
    _is_real_number_track: int = PrivateAttr(default=0)

    @field_validator("stat_method")
    @classmethod
    def _validate_stat_method(cls, value):
        if value is None:
            return value
        if value not in _BIN_STATS:
            raise ValueError(f"Unsupported bin statistic: {value}")
        return value

    @field_validator("data_transform", mode="before")
    @classmethod
    def _coerce_data_transform(cls, value):
        if value is None:
            return _echo
        if isinstance(value, str):
            if value in _TRANSFORMATIONS:
                return _TRANSFORMATIONS[value]
            raise ValueError(f"Transformation is not supported ({value})")
        if callable(value):
            return value
        raise ValueError(f"Transformation is not supported ({value})")

    @field_validator("convert_nan_to_num", mode="before")
    @classmethod
    def _coerce_nan_converter(cls, value):
        if value is None:
            return _echo
        if callable(value):
            return value
        raise ValueError(
            "value of convert_nan_to_num must be None or a callable object."
        )

    def _get(self, chromosome, start, end):
        raise NotImplementedError

    @staticmethod
    def _echo(data):
        return _echo(data)

    def reset_min_val(self):
        """Remove constraints for min value."""
        self.min_val = None

    def reset_max_val(self):
        """Remove constraints for max value."""
        self.max_val = None

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(chromosome, start, end, ax, index=index, **kwargs)
        if index != 0:
            self._ax.set_xticks([])
            self._ax.set_xticks([], minor=True)
        else:
            ax.xaxis.tick_top()
        self._ax.margins(0)

    def _get_scale(self, a=1):
        def forward(x):
            x = (x >= 0) * x + (x < 0) * x * a
            return x

        def inverse(x):
            x = (x >= 0) * x + (x < 0) * x / a
            return x

        return forward, inverse

    def _post_plot_hook(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._post_plot_hook(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        y_start, y_end = self._ax.get_ylim()
        if self.min_val is not None:
            y_start = self.min_val
        if self.max_val is not None:
            y_end = self.max_val

        if self._is_real_number_track:
            if y_start < 0:
                if self.equal_space_for_pos_neg_ranges:
                    if y_end > 0:
                        forward, inverse = self._get_scale(np.abs(y_end / y_start))
                        self._yscale_func = (forward, inverse)
                        ax.set_yscale("function", functions=(forward, inverse))
                    ranges = [y_start, 0, y_end]
                else:
                    ranges = [y_start, 0, y_end]
            else:
                ranges = [y_start, y_end]
        else:
            ranges = [y_start, y_end]
        self._ax.set_ylim((y_start, y_end))
        if self.show_range:
            self._ax.yaxis.set_ticks(ranges)
        else:
            self._ax.yaxis.set_ticks([])

        if self.y_tick_format is not None:
            ticks = self._ax.get_yticks()
            new_labels = [
                self.y_tick_format.format(label) for label in self._ax.get_yticks()
            ]
            self._ax.set_yticks(ticks)
            self._ax.set_yticklabels(new_labels)

        if self.inward_yticks:
            y_ticks = self._ax.get_yticklabels()
            if len(y_ticks) > 1:
                y_ticks[0].set_verticalalignment("bottom")
                y_ticks[-1].set_verticalalignment("top")

        distance_cutoff = max(0.01 * end - start, 1)
        if self.max_val is not None or self.min_val is not None:
            for line in self._ax.get_lines():
                t = line.get_xydata()
                x = t[:, 0]
                y = t[:, 1]
                if self.max_val is not None:
                    to_be_masked = np.logical_and(y > self.max_val, y > y_end)
                    n_to_be_masked = to_be_masked.sum()
                    if n_to_be_masked > 0:
                        self._ax.scatter(
                            x[to_be_masked],
                            [y_end] * n_to_be_masked,
                            marker="_",
                            color="black",
                            s=2,
                            zorder=100,
                        )
                    if self.label_masked_peak:
                        from scipy.signal import find_peaks

                        if distance_cutoff >= 1:
                            peaks, _ = find_peaks(
                                y,
                                rel_height=1,
                                height=self.max_val,
                                distance=distance_cutoff,
                            )
                            texts = []
                            ha_choices = ("right", "left")
                            for i, _x in enumerate(peaks):
                                X = x[_x]
                                Y = y[_x]
                                if self.overflow_label_format is not None:
                                    s = self.overflow_label_format.format(Y)
                                else:
                                    s = "{:.2f}".format(Y)
                                if not self.overflow_label_auto_adjust:
                                    texts.append(
                                        self._ax.text(
                                            X,
                                            self.max_val,
                                            s,
                                            va="bottom",
                                            ha=ha_choices[i % 2],
                                        )
                                    )
                                else:
                                    texts.append(
                                        self._ax.text(X, self.max_val, s, va="bottom")
                                    )
                            if self.overflow_label_auto_adjust and len(texts) > 0:
                                try:
                                    from adjustText import adjust_text

                                    adjust_text(texts)
                                except ImportError:
                                    pass

                if self.min_val is not None:
                    to_be_masked = np.logical_and(y < self.min_val, y < y_start)
                    n_to_be_masked = to_be_masked.sum()
                    if n_to_be_masked > 0:
                        self._ax.scatter(
                            x[to_be_masked],
                            [y_start] * n_to_be_masked,
                            marker="_",
                            color="black",
                            s=2,
                            zorder=100,
                        )

                    if self.label_masked_peak:
                        from scipy.signal import find_peaks

                        peaks, _ = find_peaks(
                            -1 * y,
                            rel_height=1,
                            height=-1 * self.min_val,
                            distance=distance_cutoff,
                        )
                        texts = []
                        ha_choices = ("right", "left")
                        for i, _x in enumerate(peaks):
                            X = x[_x]
                            Y = y[_x]
                            if self.overflow_label_format is not None:
                                s = self.overflow_label_format.format(Y)
                            else:
                                s = "{:.2f}".format(Y)
                            if not self.overflow_label_auto_adjust:
                                texts.append(
                                    self._ax.text(
                                        X,
                                        self.min_val,
                                        s,
                                        va="bottom",
                                        ha=ha_choices[i % 2],
                                    )
                                )
                            else:
                                texts.append(
                                    self._ax.text(X, self.min_val, s, va="bottom")
                                )
                        if self.overflow_label_auto_adjust and len(texts) > 0:
                            try:
                                from adjustText import adjust_text

                                adjust_text(texts)
                            except Exception:
                                pass

    def _merge_redundant_values(self, x: np.ndarray, y: np.ndarray) -> list:
        keep_idx = [0]
        for i in range(1, x.shape[0] - 1):
            if y[i + 1] == y[i] and y[i] == y[i - 1]:
                continue
            else:
                keep_idx.append(i)
        keep_idx.append(x.shape[0] - 1)
        return keep_idx


class DynamicValueTrack(NumericalTrack):
    """Show numerical values assigned in code rather than loaded from a file."""

    track: str = Field(default="", description="Placeholder")

    _values: Any = PrivateAttr(default=None)

    def __init__(self, track: str = "", **data: Any) -> None:
        super().__init__(track=track, **data)

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, value):
        self._values = value

    def _get(self, chromosome, start, end):
        xvalues = np.arange(start, end, step=1)
        if len(xvalues) != len(self.values):
            raise ValueError(
                "The length of the region (end-start) is different from values' length."
            )
        return xvalues, self.values

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        x, y = self._get(chromosome=chromosome, start=start, end=end)
        self._ax.plot(
            x, y, color=self.color, linewidth=self.line_width, alpha=self.alpha
        )
        self._ax.fill_between(
            x, y, 0, facecolor=self.color, alpha=self.alpha, lw=self.line_width
        )


class DualAxisTrack(Track):
    """Compose two existing tracks into one subplot with dual y-axes."""

    left_track: Track = Field(description="Track drawn on the left y-axis")
    right_track: Track = Field(description="Track drawn on the right y-axis")

    _right_ax: Any = PrivateAttr(default=None)

    def __init__(self, left_track: Track, right_track: Track, **data: Any) -> None:
        super().__init__(left_track=left_track, right_track=right_track, **data)

    def layout_height(self) -> float:
        return max(self.left_track.layout_height(), self.right_track.layout_height())

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        inward_ticks = kwargs.get("inward_ticks", False)
        self.left_track._pre_plot_hook(
            chromosome=chromosome,
            start=start,
            end=end,
            inward_ticks=inward_ticks,
        )
        self.right_track._pre_plot_hook(
            chromosome=chromosome,
            start=start,
            end=end,
            inward_ticks=inward_ticks,
        )

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        self.left_track._draw_track(
            chromosome=chromosome,
            start=start,
            end=end,
            ax=ax,
            index=index,
            **kwargs,
        )
        self._ax = self.left_track._ax

        self._right_ax = ax.twinx()
        right_kwargs = dict(kwargs)
        right_kwargs["hide_coordinates"] = True
        self.right_track._draw_track(
            chromosome=chromosome,
            start=start,
            end=end,
            ax=self._right_ax,
            index=1,
            **right_kwargs,
        )
        self._right_ax.set_xlim(self._ax.get_xlim())
        self._right_ax.yaxis.set_label_position("right")
        self._right_ax.yaxis.tick_right()
        self._right_ax.spines["left"].set_visible(False)

    def _post_plot_hook(self, chromosome, start, end, ax, index=1, **kwargs):
        self.left_track._post_plot_hook(
            chromosome=chromosome,
            start=start,
            end=end,
            ax=self.left_track._ax,
            index=index,
            **kwargs,
        )
        self.right_track._post_plot_hook(
            chromosome=chromosome,
            start=start,
            end=end,
            ax=self.right_track._ax,
            index=1,
            **kwargs,
        )
        self._right_ax.set_xlim(self.left_track._ax.get_xlim())

    def set_highlight_regions(self, starts, ends, colors=(), alpha_vals=()):
        self.left_track.set_highlight_regions(starts, ends, colors, alpha_vals)
        self.right_track.set_highlight_regions(starts, ends, colors, alpha_vals)

    def add_highlight_region(self, start, end):
        self.left_track.add_highlight_region(start, end)
        self.right_track.add_highlight_region(start, end)

    def remove_highlight(self):
        self.left_track.remove_highlight()
        self.right_track.remove_highlight()
