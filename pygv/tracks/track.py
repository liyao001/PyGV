from abc import abstractmethod
from typing import Any, Callable, ClassVar, Dict, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from pydantic import BaseModel, ConfigDict, Field, field_validator

from pygv.errors.DataIntegrity import InvaildRegion
from pygv.errors.Implementation import UnimplementedBinStat, UnimplementedTransformation


class _ConfigModel(BaseModel):
    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        extra="forbid",
        populate_by_name=True,
        validate_assignment=True,
    )

    def __init__(self, **data: Any):
        super().__init__(**self._populate_config_defaults(data))

    @classmethod
    def _populate_config_defaults(cls, data):
        data = dict(data)
        defaults = {}
        for base in reversed(cls.mro()):
            if (
                isinstance(base, type)
                and issubclass(base, _ConfigModel)
                and base is not _ConfigModel
                and base.__name__.endswith("Config")
            ):
                for name, field in base.model_fields.items():
                    if field.is_required():
                        continue
                    alias = field.alias
                    if name in data or (alias is not None and alias in data):
                        continue
                    defaults[name] = field.get_default(call_default_factory=True)
        defaults.update(data)
        return defaults


class TrackConfig(_ConfigModel):
    name: str = Field(default="", description="Track label shown on the y-axis.")
    line_width: float = Field(
        default=1, ge=0, description="Default line width for rendered track elements."
    )
    height: float = Field(
        default=1, gt=0, description="Relative height of the track in the viewer layout."
    )
    color: Optional[Any] = Field(
        default="#A1A1A1", description="Default Matplotlib color for track elements."
    )
    edge_color: Optional[Any] = Field(
        default="#6E6E6E", description="Default Matplotlib edge color for track elements."
    )
    font_color: Optional[Any] = Field(
        default="black", description="Matplotlib color used for track text."
    )
    alpha: float = Field(
        default=0.8,
        ge=0,
        le=1,
        description="Opacity for patches and filled visual elements.",
    )
    font_size: float = Field(
        default_factory=lambda: plt.rcParams["font.size"],
        gt=0,
        description="Font size used for track labels and text.",
    )
    y_tick_format: Optional[str] = Field(
        default=None, description="Format string used for y-axis tick labels."
    )
    y_label_rotation: Union[float, str] = Field(
        default="horizontal",
        description="Rotation of the y-axis label: 'horizontal', 'vertical', or a numeric angle.",
    )
    y_label_ha: str = Field(
        default="right", description="Horizontal alignment for the y-axis label."
    )
    y_label_va: str = Field(
        default="center", description="Vertical alignment for the y-axis label."
    )
    inward_yticks: bool = Field(
        default=False, description="Adjust end y-tick labels inward to reduce overlap."
    )

    @field_validator("color", "edge_color", "font_color")
    def _validate_color(cls, value):
        if value is None or is_color_like(value):
            return value
        raise ValueError(f"Invalid color value: {value}")

    @field_validator("y_label_rotation")
    def _validate_y_label_rotation(cls, value):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return value
        if value in {"vertical", "horizontal"}:
            return value
        raise ValueError(
            "y_label_rotation must be 'vertical', 'horizontal', or a numeric angle."
        )


class AnnotationTrackConfig(TrackConfig):
    patch_height: float = Field(
        default=1, gt=0, description="Height of annotation patches such as exons or blocks."
    )
    allowed_feature_lanes: Optional[int] = Field(
        default=None,
        gt=0,
        description="Maximum number of annotation lanes to draw; None allows all lanes.",
    )
    font_box_alpha: float = Field(
        default=0.75,
        ge=0,
        le=1,
        description="Opacity of text background boxes used for annotation labels.",
    )
    lane_space: float = Field(
        default=0.25, ge=0, description="Extra vertical spacing between annotation lanes."
    )
    features_per_lane: int = Field(
        default=3, gt=0, description="Maximum number of features grouped into each lane."
    )
    line_color: Optional[Any] = Field(
        default="black", description="Matplotlib color used for annotation connector lines."
    )
    arrow_interval: float = Field(
        default=5, gt=0, description="Spacing between strand direction arrows."
    )
    padding_left: float = Field(
        default=0, ge=0, description="Extra left-side spacing used when placing features."
    )
    padding_right: float = Field(
        default=0, ge=0, description="Extra right-side spacing used when placing features."
    )
    show_name: bool = Field(default=True, description="Whether to draw feature names.")
    hide_visual_dup: bool = Field(
        default=False,
        description="Whether to hide visually duplicated features in the current view.",
    )

    @field_validator("line_color")
    def _validate_line_color(cls, value):
        if value is None or is_color_like(value):
            return value
        raise ValueError(f"Invalid color value: {value}")


class NumericalTrackConfig(TrackConfig):
    min_val: Optional[float] = Field(
        default=None, description="Lower y-axis limit for numerical signals."
    )
    max_val: Optional[float] = Field(
        default=None, description="Upper y-axis limit for numerical signals."
    )
    show_range: bool = Field(default=True, description="Whether to show y-axis range ticks.")
    n_bins: Optional[int] = Field(
        default=None, gt=0, description="Number of bins used to summarize numerical values."
    )
    stat_method: Optional[str] = Field(
        default=None, description="Statistic used when summarizing values into bins."
    )
    data_transform: Optional[Union[str, Callable]] = Field(
        default=None,
        alias="transformation",
        description="Transformation applied to numerical values before plotting.",
    )
    convert_nan_to_num: Optional[Callable] = Field(
        default=np.nan_to_num,
        description="Callable used to convert NaN values before plotting.",
    )
    scale: float = Field(
        default=1, description="Multiplicative scaling factor applied to numerical values."
    )
    label_masked_peak: bool = Field(
        default=True, description="Whether to label peaks clipped by min_val or max_val."
    )
    overflow_label_format: Optional[str] = Field(
        default="{:.1f}", description="Format string used for clipped peak labels."
    )
    overflow_label_auto_adjust: bool = Field(
        default=False,
        description="Whether to automatically adjust clipped peak label placement.",
    )

    @field_validator("stat_method")
    def _validate_stat_method(cls, value):
        if value is None:
            return value
        np_supported_methods = {"mean", "std", "median", "count", "sum", "min", "max"}
        if value in np_supported_methods:
            return value
        raise ValueError(f"Unsupported stat_method: {value}")

    @field_validator("convert_nan_to_num")
    def _validate_nan_converter(cls, value):
        if value is None or callable(value):
            return value
        raise ValueError("convert_nan_to_num must be None or a callable object.")

    @field_validator("data_transform")
    def _validate_data_transform(cls, value):
        if value is None or callable(value):
            return value
        supported_transformations = {
            "ln",
            "asinh",
            "log2",
            "log10",
            "log1p",
            "rln",
            "rlog2",
            "rlog10",
            "rlog1p",
        }
        if value in supported_transformations:
            return value
        raise ValueError(str(UnimplementedTransformation(value)))


class DynamicValueTrackConfig(NumericalTrackConfig):
    values: Any = Field(default=None, description="Dynamic values to plot for the track.")


class Track(TrackConfig):
    """
    Generic Track

    Parameters
    ----------
    kwargs : dict
        name : str
            Name of the track
        line_width : numeric
            The default width for lines
        height : int
            Height of the track (unit, relative measurement)
        color : color_like
            Default color, #A1A1A1
        edge_color : color_like
            Edge color, #6E6E6E
        font_color : color_like
            Font color, black
        alpha : float
            Alpha of patches
        font_size : float
            Font size
        y_tick_format : str
            String format for ticks on y-axis. For example: `{:.1f}` (only keep one digit)
        highlight_start : int
            Start loc of highlight region
        highlight_end : int
            End loc of highlight region
        highlight_color : color_like
            Highlight color
        highlight_alpha : float
            Alpha for highlighting
        y_label_rotation : str or float
            Rotation of y-axis' label, by default, vertical.
        y_label_ha : str
            Horizontal alignment about label for y-axis
    """

    _FIELD_PRIVATE_ATTRS: ClassVar[Dict[str, str]] = {
        "name": "_name",
        "line_width": "_line_width",
        "height": "_height",
        "color": "_color",
        "edge_color": "_edge_color",
        "font_color": "_font_color",
        "alpha": "_alpha",
        "font_size": "_font_size",
        "y_tick_format": "_y_tick_format",
        "y_label_rotation": "_y_label_rotation",
        "y_label_ha": "_y_label_ha",
        "y_label_va": "_y_label_va",
        "inward_yticks": "_inward_yticks",
    }

    def __getattribute__(self, name):
        private_attrs = object.__getattribute__(self, "_FIELD_PRIVATE_ATTRS")
        private_name = private_attrs.get(name)
        if private_name is not None:
            values = object.__getattribute__(self, "__dict__")
            if private_name in values:
                return values[private_name]
        return super().__getattribute__(name)

    def __setattr__(self, name, value):
        if name in type(self).model_fields:
            super().__setattr__(name, value)
            self._sync_private_field(name)
        else:
            object.__setattr__(self, name, value)

    def _sync_private_field(self, name):
        private_name = self._FIELD_PRIVATE_ATTRS.get(name)
        if private_name is not None and name in self.__dict__:
            object.__setattr__(self, private_name, self.__dict__[name])

    def _sync_private_fields(self):
        for name in self._FIELD_PRIVATE_ATTRS:
            self._sync_private_field(name)

    def dict(self, *args, **kwargs):
        if kwargs.get("include") is None:
            kwargs["include"] = set(type(self).model_fields)
        return super().model_dump(*args, **kwargs)

    def model_dump(self, *args, **kwargs):
        if kwargs.get("include") is None:
            kwargs["include"] = set(type(self).model_fields)
        return super().model_dump(*args, **kwargs)

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self._sync_private_fields()
        self._ax = None

        # highlight spans
        self._highlight_starts = []
        self._highlight_ends = []
        self._highlight_colors = []
        self._highlight_alphas = []

    @property
    def alpha(self):
        """
        Alpha of patches
        """
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        if 0 <= value <= 1:
            self._alpha = value
        else:
            raise ValueError("alpha must be between 0 and 1")

    @property
    def name(self):
        """
        Name of the track
        """
        return self._name

    @name.setter
    def name(self, value):
        try:
            self._name = str(value)
        except:
            pass

    @property
    def color(self):
        """
        Default color
        """
        return self._color

    @color.setter
    def color(self, value):
        try:
            if is_color_like(value) or value is None:
                self._color = str(value)
            else:
                print(f"Invalid color value: {value}")
        except:
            pass

    @property
    def edge_color(self):
        """
        Edge color
        """
        return self._edge_color

    @edge_color.setter
    def edge_color(self, value):
        try:
            if is_color_like(value) or value is None:
                self._edge_color = str(value)
            else:
                print(f"Invalid color value: {value}")
        except:
            pass

    @property
    def font_color(self):
        """
        Font color
        """
        return self._font_color

    @font_color.setter
    def font_color(self, value):
        try:
            if is_color_like(value) or value is None:
                self._font_color = value
            else:
                print(f"Invalid color value {value}")
        except:
            pass

    @property
    def font_size(self):
        """
        Font size
        """
        return self._font_size

    @font_size.setter
    def font_size(self, value):
        try:
            self._font_size = float(value)
        except:
            pass

    @property
    def line_width(self):
        """
        Line width
        """
        return self._line_width

    @line_width.setter
    def line_width(self, value):
        try:
            self._line_width = float(value)
        except:
            pass

    @property
    def height(self):
        """
        Height of the track (unit, relative measurement)
        """
        return self._height

    @height.setter
    def height(self, value):
        if value > 0:
            self._height = value
        else:
            raise ValueError("height must be greater than 0")

    @property
    def y_tick_format(self):
        """
        String format for ticks on y-axis
        """
        return self._y_tick_format

    @y_tick_format.setter
    def y_tick_format(self, value):
        self._y_tick_format = value

    @property
    def y_label_rotation(self):
        """
        Rotation of y-axis' label, by default, vertical.
        """
        return self._y_label_rotation

    @y_label_rotation.setter
    def y_label_rotation(self, value):
        self._y_label_rotation = value

    @property
    def y_label_ha(self):
        """
        Set the horizontal alignment
        """
        return self._y_label_ha

    @y_label_ha.setter
    def y_label_ha(self, value):
        self._y_label_ha = value

    @property
    def y_label_va(self):
        """
        Set the vertical alignment
        """
        return self._y_label_va

    @y_label_va.setter
    def y_label_va(self, value):
        self._y_label_va = value

    @property
    def inward_yticks(self):
        """
        Plot y-ticks strictly inside each track.
        If you want to apply inward_yticks to all tracks,
        you can set `inward_yticks=True` when creating the `GenomeViewer`,
        like `GenomeViewer(inward_yticks=True)`

        Examples
        --------

        .. plot:: ../examples/plot_inward_yticks.py
        """
        return self._inward_yticks

    @inward_yticks.setter
    def inward_yticks(self, value):
        if value is not None:
            self._inward_yticks = bool(value)

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        """
        This method will be called before calling the :func:`~pygv.tracks.track.Track.draw_track` method.
        For now, it sets the default font size

        Parameters
        ----------
        chromosome : str
            chromosome
        start : int
            start of the ROI, 0-based
        end : int
            end of the ROI, 0-based
        kwargs

        Returns
        -------

        """
        if self._font_size is None:
            self._font_size = plt.rcParams["font.size"]
        inward_ticks = kwargs.pop("inward_ticks", False)
        if inward_ticks is not None:
            self.inward_yticks = inward_ticks

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw track

        Parameters
        ----------
        chromosome : str
            chromosome
        start : int
            start of the ROI, 0-based
        end : int
            end of the ROI, 0-based
        ax : matplotlib.axes.Axes
            matplotlib.axes.Axes to plot
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        self._ax = ax
        n_ticks = kwargs.get("n_ticks", None)
        hide_coords = kwargs.get("hide_coordinates", False)
        hide_chr_name = kwargs.get("hide_chromosome_name", None)
        if start <= end:
            self._ax.set_xlim((start, end))
        else:
            raise InvaildRegion("start of the region must be smaller than the end")
        if index != 0 or hide_coords:
            self._ax.spines["top"].set_visible(False)
            # remove major ticks
            self._ax.set_xticks([])
            # remove minor ticks
            self._ax.set_xticks([], minor=True)
        else:
            # plot coordinates
            self._ax.xaxis.set_ticks_position("top")
            self._ax.spines["top"].set_position(("outward", 10))
            self._ax.spines["top"].set_linewidth(2)
            if n_ticks is not None:
                ticks = np.linspace(start, end, n_ticks, dtype=int)
            else:
                ticks = [t for t in self._ax.get_xticks() if start <= t <= end]

            if ticks[-1] - ticks[1] <= 1e3:
                labels = [f"{x:,.0f}"
                          for x in ticks]
                labels[-1] += " bp"

            elif ticks[-1] - ticks[1] <= 4e5:
                labels = [f"{x / 1000.0:,.0f}"
                          for x in ticks]
                labels[-1] += " Kb"

            else:
                labels = [f"{x / 1000000.0:,.1f} "
                          for x in ticks]
                labels[-1] += " Mbp"

            if not hide_chr_name:
                self._ax.set_title(chromosome)
            self._ax.set_xticks(ticks)
            self._ax.set_xticklabels(labels)

        self._ax.spines["bottom"].set_visible(False)
        self._ax.spines["right"].set_visible(False)

        if self.name is not None:
            self._ax.set_ylabel(
                self.name,
                rotation=self.y_label_rotation,
                ha=self.y_label_ha,
                va=self.y_label_va,
            )

    def set_highlight_regions(self, starts, ends, colors=(), alpha_vals=()):
        """
        Set highlight region

        Parameters
        ----------
        starts : list of numeric values
            Start positions of the highlight zones
        ends : list of numeric values
            End positions of the highlight zones
        colors : tuple
            Leave it as an empty tuple if you want to use the default color.
            If you only give one color, it will be applied to all regions; otherwise, you should specify
            colors for each region.
        alpha_vals : tuple
            Leave it as an empty tuple if you want to use the default transparency level (0.5).
            If you only give one value, it will be applied to all regions; otherwise, you should specify
            transparency values for each region.

        Returns
        -------

        """
        try:
            if isinstance(starts, int) or isinstance(starts, float):
                starts = [
                    starts,
                ]
            if isinstance(ends, int) or isinstance(ends, float):
                ends = [
                    ends,
                ]
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
        except:
            pass

    def add_highlight_region(self, start, end):
        """
        Set highlight region

        Parameters
        ----------
        start : numeric
            Start position of the highlight zone
        end : numeric
            End position of the highlight zone

        Returns
        -------

        """
        try:
            start = float(start)
            end = float(end)
            self._highlight_starts.append(start)
            self._highlight_ends.append(end)
        except:
            pass

    def remove_highlight(self):
        """
        Remove highlight zone

        Returns
        -------

        """
        self._highlight_starts = []
        self._highlight_ends = []

    def _post_plot_hook(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        The hook is called after :func:`~pygv.tracks.track.Track.draw_track` method,
        the method here checks highlight settings and if there's any,
        they will be highlighted with highlight_color and highlight_alpha

        Parameters
        ----------
        chromosome : str
            chromosome
        start : int
            start of the ROI, 0-based
        end : int
            end of the ROI, 0-based
        ax : matplotlib.axes.Axes
            matplotlib.axes.Axes to plot
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
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


class AnnotationTrack(Track, AnnotationTrackConfig):
    """
    Annotation track

    Parameters
    ----------
    track : str
    kwargs : dict
        allowed_feature_lanes : None or int
            :attr:`allowed_feature_lanes`
        height : float
            :attr:`height`
        arrow_interval :
            :attr:`arrow_interval`
        features_per_lane : int
            :attr:`features_per_lane`
        font_box_alpha : float
            :attr:`font_box_alpha`
        lane_space : float
            :attr:`lane_space`
        line_color :
            :attr:`line_color`
        padding_left : int
            :attr:`padding_left`
        padding_right : int
            :attr:`padding_right`
        patch_height : int
            :attr:`patch_height`
        hide_visual_dup : bool
            :attr:`hide_visual_dup`

    """

    _FIELD_PRIVATE_ATTRS: ClassVar[Dict[str, str]] = {
        **Track._FIELD_PRIVATE_ATTRS,
        "patch_height": "_patch_height",
        "allowed_feature_lanes": "_allowed_feature_lanes",
        "font_box_alpha": "_font_box_alpha",
        "lane_space": "_lane_space",
        "features_per_lane": "_features_per_lane",
        "line_color": "_line_color",
        "arrow_interval": "_arrow_interval",
        "padding_left": "_padding_left",
        "padding_right": "_padding_right",
        "show_name": "_show_name",
        "hide_visual_dup": "_hide_visual_dup",
    }

    def __getattribute__(self, name):
        if name == "height":
            values = object.__getattribute__(self, "__dict__")
            if "_lane_registries" in values and "_height" in values:
                return max(len(values["_lane_registries"]), 1) * values["_height"]
        return super().__getattribute__(name)

    @property
    def patch_height(self):
        """
        Height of patches (for exons/blocks)
        """
        return self._patch_height

    @patch_height.setter
    def patch_height(self, value):
        try:
            self._patch_height = float(value)
        except:
            pass

    @property
    def allowed_feature_lanes(self):
        """
        Max amount of feature lanes to be plotted. For example, if a region has 12 overlapping features,
        to make sure all features can be clearly rendered, these features will be plotted into 12 separate lanes.
        If you set the value to be smaller than 12 (say `2`), then you will only see two lanes in the end.

        Examples
        --------

        .. plot:: ../examples/plot_allowed_feature_lanes.py
        """
        return self._allowed_feature_lanes

    @allowed_feature_lanes.setter
    def allowed_feature_lanes(self, value):
        try:
            if value is None:
                self._allowed_feature_lanes = None
            else:
                self._allowed_feature_lanes = int(value)
        except:
            pass

    @property
    def font_box_alpha(self):
        """
        Transparent/alpha for text boxes labeling gene names
        """
        return self._font_box_alpha

    @font_box_alpha.setter
    def font_box_alpha(self, value):
        try:
            self._font_box_alpha = float(value)
        except:
            pass

    @property
    def lane_space(self):
        """
        Extra spaces between lanes
        """
        return self._lane_space

    @lane_space.setter
    def lane_space(self, value):
        try:
            self._lane_space = float(value)
        except:
            pass

    @property
    def features_per_lane(self):
        """
        Features per lane
        """
        return self._features_per_lane

    @features_per_lane.setter
    def features_per_lane(self, value):
        try:
            self._features_per_lane = int(value)
        except:
            pass

    @property
    def line_color(self):
        """
        Line color
        """
        return self._line_color

    @line_color.setter
    def line_color(self, value):
        try:
            self._line_color = value
        except:
            pass

    @property
    def arrow_interval(self):
        """
        Intervals between arrows
        """
        return self._arrow_interval

    @arrow_interval.setter
    def arrow_interval(self, value):
        try:
            self._arrow_interval = float(value)
        except:
            pass

    @property
    def padding_left(self):
        """
        To ensure that feature names do not overlap with one another, you can introduce additional
        padding spaces on the left side of each feature. When setting an integer value (let's call
        it :math:`x`) for this property, features will be placed in separate lanes if the distance between
        them is less than x. Alternatively, if you opt for a float value between 0 and 1 (designated as :math:`f`),
        the required spacing will be a fraction of the current visible region's length (denoted as :math:`l`),
        making the final spacing requirement equal to :math:`l\\times f`.
        """
        return self._padding_left

    @padding_left.setter
    def padding_left(self, value):
        try:
            self._padding_left = float(value)
        except:
            pass

    @property
    def show_name(self):
        """
        By default, PyGV prints the names of genomic regions if available.
        This behavior can be changed by assigning :code:`False` to this property.
        """
        return self._show_name

    @show_name.setter
    def show_name(self, value):
        try:
            self._show_name = bool(value)
        except:
            pass

    @property
    def hide_visual_dup(self):
        """
        Hide features which are "duplicates" to other features in current window (only one will be kept)
        """
        return self._hide_visual_dup

    @hide_visual_dup.setter
    def hide_visual_dup(self, value):
        try:
            self._hide_visual_dup = bool(value)
        except:
            pass

    @property
    def height(self):
        return max(len(self._lane_registries), 1) * self._height

    @height.setter
    def height(self, value):
        self._height = value

    def __init__(self, track, **kwargs):
        super(AnnotationTrack, self).__init__(**kwargs)

        self._plot_thickness = 0
        self._plot_block = 0
        self._small_relative = 0

        # override defaults
        if self.color is None:
            self.color = "#A1A1A1"
        if self.edge_color is None:
            self.edge_color = "#6E6E6E"

        # lane manager
        self._lane_registries = []

    def _plot_gene_direction(self, ax, xpos, ypos, strand, **kwargs):
        """
        Draws a broken line with 2 parts:
        For strand = +:  > For strand = -: <
        :param xpos:
        :param ypos:
        :param strand:
        :
        :return: None
        """
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
            Line2D(xdata, ydata, color=self._line_color, linewidth=self.line_width)
        )


class NumericalTrack(Track, NumericalTrackConfig):
    """
    Numerical track

    Parameters
    ----------
    kwargs : dict
        max_val : int, optional
            Maximum value to be plotted. By default, all signals are plotted.
        min_val : int, optional
            Minimum value to be plotted. By default, all signals are plotted.
        label_masked_peak : bool
            Whether or not to labelled capped signals.
        overflow_label_format : str
            String format for labeling overflow loci
        overflow_label_auto_adjust : bool
            Switch controlling the automatic placement of text labels for overflow signals
    """

    _FIELD_PRIVATE_ATTRS: ClassVar[Dict[str, str]] = {
        **Track._FIELD_PRIVATE_ATTRS,
        "min_val": "_min_val",
        "max_val": "_max_val",
        "show_range": "_show_range",
        "n_bins": "_n_bins",
        "stat_method": "_stat_method",
        "data_transform": "_data_transform",
        "convert_nan_to_num": "_convert_nan_to_num",
        "scale": "_scale",
        "label_masked_peak": "_label_masked_peak",
        "overflow_label_format": "_overflow_label_format",
        "overflow_label_auto_adjust": "_overflow_label_auto_adjust",
    }

    def _sync_private_field(self, name):
        if name == "data_transform":
            object.__setattr__(
                self,
                "_data_transform",
                self._resolve_data_transform(self.__dict__[name]),
            )
            return
        if name == "convert_nan_to_num":
            value = self.__dict__[name]
            object.__setattr__(
                self,
                "_convert_nan_to_num",
                self._echo if value is None else value,
            )
            return
        super()._sync_private_field(name)

    @staticmethod
    def _resolve_data_transform(value):
        if value is None:
            return NumericalTrack._echo
        if callable(value):
            return value
        transformations = {
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
        try:
            return transformations[value]
        except KeyError:
            raise UnimplementedTransformation(value)

    @abstractmethod
    def _get(self, chromosome, start, end):
        pass

    @staticmethod
    def _echo(data):
        return data

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super(NumericalTrack, self)._draw_track(
            chromosome, start, end, ax, index=index, **kwargs
        )
        if index != 0:
            # remove major ticks
            self._ax.set_xticks([])
            # remove minor ticks
            self._ax.set_xticks([], minor=True)
        else:
            ax.xaxis.tick_top()
        self._ax.margins(0)

    def _get_scale(self, a=1):
        """
        Source: https://stackoverflow.com/questions/53699677/matplotlib-different-scale-on-negative-side-of-the-axis

        Parameters
        ----------
        a

        Returns
        -------

        """

        def forward(x):
            x = (x >= 0) * x + (x < 0) * x * a
            return x

        def inverse(x):
            x = (x >= 0) * x + (x < 0) * x / a
            return x

        return forward, inverse

    def _post_plot_hook(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Post-plot hooks

        Parameters
        ----------
        chromosome : str
            chromosome
        start : int

        end : int

        ax : matplotlib.axes.Axes

        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs

        Returns
        -------

        """
        super(NumericalTrack, self)._post_plot_hook(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        y_start, y_end = self._ax.get_ylim()
        if self._min_val is not None:
            y_start = self._min_val
        if self._max_val is not None:
            y_end = self._max_val

        if self.is_real_number_track:
            if y_start < 0:
                if self.equal_space_for_pos_neg_ranges:
                    # only apply scale adjustment when there are both positive and negative data values
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

        if self._y_tick_format is not None:
            ticks = self._ax.get_yticks()
            new_labels = [
                self._y_tick_format.format(label) for label in self._ax.get_yticks()
            ]
            self._ax.set_yticks(ticks)
            self._ax.set_yticklabels(new_labels)

        if self.inward_yticks:
            # adjust the placements of the first and last ytick to avoid overlap
            y_ticks = self._ax.get_yticklabels()
            if len(y_ticks) > 1:
                y_ticks[0].set_verticalalignment("bottom")
                y_ticks[-1].set_verticalalignment("top")

        # add bars to show overflowed signals
        distance_cutoff = max(0.01 * end - start, 1)
        if self._max_val is not None or self._min_val is not None:
            for line in self._ax.get_lines():
                t = line.get_xydata()
                x = t[:, 0]
                y = t[:, 1]
                if self._max_val is not None:
                    to_be_masked = np.logical_and(y > self._max_val, y > y_end)
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

                        # find peaks in signal tracks which are higher than the threshold
                        # then add annotation to show their values
                        # require the distance between peaks to be away from their neighbours
                        # for at least 1% of the overall window
                        if distance_cutoff >= 1:
                            peaks, _ = find_peaks(
                                y,
                                rel_height=1,
                                height=self._max_val,
                                distance=distance_cutoff,
                            )
                            texts = []
                            ha_choices = ("right", "left")
                            for i, _x in enumerate(peaks):
                                X = x[_x]
                                Y = y[_x]
                                if self._overflow_label_format is not None:
                                    s = self._overflow_label_format.format(Y)
                                else:
                                    s = "{:.2f}".format(Y)
                                if not self._overflow_label_auto_adjust:
                                    texts.append(
                                        self._ax.text(
                                            X,
                                            self._max_val,
                                            s,
                                            va="bottom",
                                            ha=ha_choices[i % 2],
                                        )
                                    )
                                else:
                                    texts.append(
                                        self._ax.text(X, self._max_val, s, va="bottom")
                                    )
                            if self._overflow_label_auto_adjust and len(texts) > 0:
                                try:
                                    from adjustText import adjust_text

                                    adjust_text(texts)
                                except ImportError:
                                    pass

                if self._min_val is not None:
                    to_be_masked = np.logical_and(y < self._min_val, y < y_start)
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

                        # find peaks in signal tracks which are higher than the threshold
                        # then add annotation to show their values
                        # require the distance between peaks to be away from their neighbours
                        # for at least 1% of the overall window
                        peaks, _ = find_peaks(
                            -1 * y,
                            rel_height=1,
                            height=-1 * self._min_val,
                            distance=distance_cutoff,
                        )
                        texts = []
                        ha_choices = ("right", "left")
                        for i, _x in enumerate(peaks):
                            X = x[_x]
                            Y = y[_x]
                            if self._overflow_label_format is not None:
                                s = self._overflow_label_format.format(Y)
                            else:
                                s = "{:.2f}".format(Y)
                            if not self._overflow_label_auto_adjust:
                                texts.append(
                                    self._ax.text(
                                        X,
                                        self._min_val,
                                        s,
                                        va="bottom",
                                        ha=ha_choices[i % 2],
                                    )
                                )
                            else:
                                texts.append(
                                    self._ax.text(X, self._min_val, s, va="bottom")
                                )
                        if self._overflow_label_auto_adjust and len(texts) > 0:
                            try:
                                from adjustText import adjust_text

                                adjust_text(texts)
                            except:
                                pass

    @property
    def scale(self):
        """
        Normalization factors for signals, you can set this value to normalize densities by RPM, etc.
        """
        return self._scale

    @scale.setter
    def scale(self, value):
        try:
            self._scale = float(value)
        except:
            pass

    def __init__(self, **kwargs):
        super(NumericalTrack, self).__init__(**kwargs)
        self.is_real_number_track = 0
        self._yscale_func = None

    @property
    def equal_space_for_pos_neg_ranges(self):
        """
        Set it as `True` to force data range to be independently
        """
        return self._equal_space_for_pos_neg_ranges

    @equal_space_for_pos_neg_ranges.setter
    def equal_space_for_pos_neg_ranges(self, value):
        if value:
            self._equal_space_for_pos_neg_ranges = 1
        elif value == 0 or value is False:
            self._equal_space_for_pos_neg_ranges = 0
        else:
            raise ValueError("draw_y_independently must be either 0/False or 1/True")

    @property
    def min_val(self):
        """
        Min value for the y-axis. If the signal values are smaller than `min_val`, they will be capped.
        """
        return self._min_val

    @min_val.setter
    def min_val(self, value):
        try:
            self._min_val = float(value)
        except:
            pass

    def reset_min_val(self):
        """
        Remove constraints for min value

        Returns
        -------

        """
        self._min_val = None

    @property
    def max_val(self):
        """
        Max value for the y-axis. If the signal values are greater than `min_val`, they will be capped.
        """
        return self._max_val

    @max_val.setter
    def max_val(self, value):
        try:
            self._max_val = float(value)
        except:
            pass

    @property
    def label_masked_peak(self):
        """
        If the signal values are capped, setting this value as True will write the original
        values near the cap signs.
        """
        return self._label_masked_peak

    @label_masked_peak.setter
    def label_masked_peak(self, value):
        self._label_masked_peak = bool(value)

    @property
    def overflow_label_format(self):
        """
        String format for labeling overflow loci
        """
        return self._overflow_label_format

    @overflow_label_format.setter
    def overflow_label_format(self, value):
        self._overflow_label_format = value

    @property
    def overflow_label_auto_adjust(self):
        """
        Switch controlling the automatic placement of text labels for overflow signals
        """
        return self._overflow_label_auto_adjust

    @overflow_label_auto_adjust.setter
    def overflow_label_auto_adjust(self, value):
        self._overflow_label_auto_adjust = value

    def reset_max_val(self):
        """
        Remove constraints for max value

        Returns
        -------

        """
        self._max_val = None

    @property
    def show_range(self):
        """
        Max value for the y-axis
        """
        return self._show_range

    @show_range.setter
    def show_range(self, value):
        try:
            self._show_range = bool(value)
        except:
            pass

    @property
    def convert_nan_to_num(self):
        """
        Nan conversion
        """
        return self._convert_nan_to_num

    @convert_nan_to_num.setter
    def convert_nan_to_num(self, value: Union[None, Callable]):
        """
        Convert nan values to numbers

        Parameters
        ----------
        value : Union[None, Callable]
            The function to mapping nan values. If set to None, the function will do nothing (echo).

        Returns
        -------

        """
        if value is None:
            self._convert_nan_to_num = self._echo
        elif callable(value):
            self._convert_nan_to_num = value
        else:
            raise ValueError(
                "value of convert_nan_to_num must be None or a callable object."
            )

    @property
    def n_bins(self):
        """
        Number of bins to apply, if a positive number is set, the window will be separated in bins and stat method will be applied, default `None` (raw signals)
        """
        return self._n_bins

    @n_bins.setter
    def n_bins(self, value):
        try:
            self._n_bins = int(value)
        except:
            pass

    @property
    def stat_method(self):
        """
        Statistical method for binning windows
        """
        return self._stat_method

    @stat_method.setter
    def stat_method(self, value):
        try:
            np_supported_methods = {
                "mean",
                "std",
                "median",
                "count",
                "sum",
                "min",
                "max",
            }
            if value is not None:
                if value in np_supported_methods:
                    self._stat_method = value
                else:
                    raise UnimplementedBinStat
            else:
                self._stat_method = None
        except:
            pass

    @property
    def data_transform(self):
        """
        Function for data transformation, currently supported values:

            * None: no function will be called, return raw values
            * "asinh": inverse hyperbolic sine function
            * "ln": natural logarithm function (log base e)
            * "log2": the binary logarithm function (log base 2)
            * "log10": the common logarithmic function (log base 10)
            * "log1p": the natural logarithm of one plus (ln(1+x))
            * function: a customized callable function
            Note: If you add `r` at the beginning of log functions, values will be :math:`-f(-x)`

        Examples
        --------

        .. plot:: ../examples/plot_data_transform.py
        """
        return self._data_transform

    @data_transform.setter
    def data_transform(self, value):
        try:
            if value is None:
                self._data_transform = self._echo
            elif type(value) is str:
                if value == "ln":
                    self._data_transform = np.log
                elif value == "asinh":
                    self._data_transform = np.arcsinh
                elif value == "log2":
                    self._data_transform = np.log2
                elif value == "log10":
                    self._data_transform = np.log10
                elif value == "log1p":
                    self._data_transform = np.log1p
                elif value == "rln":
                    self._data_transform = lambda x: -1 * np.log(-1 * x)
                elif value == "rlog2":
                    self._data_transform = lambda x: -1 * np.log2(-1 * x)
                elif value == "rlog10":
                    self._data_transform = lambda x: -1 * np.log10(-1 * x)
                elif value == "rlog1p":
                    self._data_transform = lambda x: -1 * np.log1p(-1 * x)
                else:
                    raise UnimplementedTransformation(value)
            elif callable(value):
                self._data_transform = value
            else:
                raise UnimplementedTransformation(value)
        except Exception as e:
            print(e)

    def _merge_redundant_values(self, x: np.ndarray, y: np.ndarray) -> list:
        keep_idx = [
            0,
        ]
        # find indices where consecutive values are not equal
        for i in range(1, x.shape[0] - 1):
            if y[i + 1] == y[i] and y[i] == y[i - 1]:
                continue
            else:
                keep_idx.append(i)
        # include the last index to ensure the final value is included
        keep_idx.append(x.shape[0] - 1)
        return keep_idx


class DynamicValueTrack(NumericalTrack, DynamicValueTrackConfig):
    """
    While other tracks load signal values from external files,
    DynamicValueTrack allows you to show the numerical values directly from your code.
    Track values should be assigned via the `values` property.

    Parameters
    ----------
    track : str
        Placeholder
    kwargs :

    Raises
    ------
    ValueError will be raised if the len the values property is not equal to the span of plotting region as defined as `end` - `start`

    Examples
    --------

    .. plot:: ../examples/plot_dyn_track.py
    """

    _FIELD_PRIVATE_ATTRS: ClassVar[Dict[str, str]] = {
        **NumericalTrack._FIELD_PRIVATE_ATTRS,
        "values": "_values",
    }

    def __init__(self, track: str = "", **kwargs):
        super(DynamicValueTrack, self).__init__(**kwargs)

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
        """
        Draw track

        Parameters
        ----------
        chromosome : str
            placeholder
        start : int
            placeholder
        end : int
            placeholder
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(DynamicValueTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        x, y = self._get(chromosome=chromosome, start=start, end=end)
        self._ax.plot(
            x, y, color=self.color, linewidth=self.line_width, alpha=self._alpha
        )
        # self.ax.bar(x, y, color=self.color, width=1)
        self._ax.fill_between(
            x, y, 0, facecolor=self.color, alpha=self._alpha, lw=self.line_width
        )
