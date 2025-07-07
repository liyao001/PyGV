from abc import abstractmethod
from typing import Any, Callable, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D

from pygv.errors.DataIntegrity import InvaildRegion
from pygv.errors.Implementation import UnimplementedBinStat, UnimplementedTransformation


class Track(object):
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

    def __init__(self, **kwargs: Any):
        # title of the track
        self._name = None
        self.name = kwargs.pop("name", "")
        # line width
        self._line_width = None
        self.line_width = kwargs.get("line_width", 1)
        # height
        self._height = None
        self.height = kwargs.pop("height", 1)
        # color
        self._color = "#A1A1A1"
        self.color = kwargs.pop("color", "#A1A1A1")
        self._edge_color = "#6E6E6E"
        self.edge_color = kwargs.pop("edge_color", "#6E6E6E")
        # alpha
        self._alpha = None
        self.alpha = kwargs.pop("alpha", 0.8)
        # font color
        self._font_color = "black"
        self.font_color = kwargs.pop("font_color", "black")
        # font size
        self._font_size = plt.rcParams["font.size"]
        self.font_size = kwargs.pop("font_size", plt.rcParams["font.size"])

        self._ax = None
        self._y_tick_format = None
        self.y_tick_format = kwargs.pop("y_tick_format", None)
        self._y_label_rotation = "horizontal"
        self.y_label_rotation = kwargs.pop("y_label_rotation", "horizontal")
        self._y_label_ha = "right"
        self.y_label_ha = kwargs.pop("y_label_ha", "right")
        self._y_label_va = "center"
        self.y_label_va = kwargs.pop("y_label_va", "center")
        self._inward_yticks = False

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
        self.inward_yticks = kwargs.pop("inward_ticks", False)

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


class AnnotationTrack(Track):
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
        # plot behaviour
        self._patch_height = 1
        self.patch_height = kwargs.pop("patch_height", 1)
        self._allowed_feature_lanes = None
        self.allowed_feature_lanes = kwargs.pop("allowed_feature_lanes", None)
        self._font_box_alpha = 0.75
        self.font_box_alpha = kwargs.pop("font_box_alpha", 0.75)
        self._lane_space = 0.25
        self.lane_space = kwargs.pop("lane_space", 0.25)
        self._features_per_lane = 3
        self.features_per_lane = kwargs.pop("features_per_lane", 3)
        self._line_color = "black"
        self.line_color = kwargs.pop("line_color", "black")
        self._arrow_interval = 5
        self.arrow_interval = kwargs.pop("arrow_interval", 5)
        self._padding_left = 0
        self.padding_left = kwargs.pop("padding_left", 0)
        self._padding_right = 0
        self.padding_right = kwargs.pop("padding_right", 0)
        self._show_name = True
        self.show_name = kwargs.pop("show_name", True)
        self._hide_visual_dup = False
        self.hide_visual_dup = kwargs.pop("hide_visual_dup", False)

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


class NumericalTrack(Track):
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

        # range
        self._min_val = None
        self.min_val = kwargs.pop("min_val", None)
        self._max_val = None
        self.max_val = kwargs.pop("max_val", None)
        self._show_range = True
        self.show_range = kwargs.pop("show_range", True)

        # stats
        self._n_bins = None
        self.n_bins = kwargs.pop("n_bins", None)
        self._stat_method = None
        self.stat_method = kwargs.pop("stat_method", None)

        # data transformation
        self._data_transform = None
        self.data_transform = kwargs.get("transformation", None)
        self._convert_nan_to_num = None
        self.convert_nan_to_num = kwargs.pop("convert_nan_to_num", np.nan_to_num)

        self._scale = None
        self.scale = kwargs.pop("scale", 1)
        self.is_real_number_track = 0

        self._label_masked_peak = None
        self.label_masked_peak = kwargs.pop("label_masked_peak", True)
        self._overflow_label_format = None
        self.overflow_label_format = kwargs.pop("overflow_label_format", "{:.1f}")
        self._overflow_label_auto_adjust = None
        self.overflow_label_auto_adjust = kwargs.pop(
            "overflow_label_auto_adjust", False
        )
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


class DynamicValueTrack(NumericalTrack):
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

    def __init__(self, track: str = "", **kwargs):
        super(DynamicValueTrack, self).__init__(**kwargs)

        self._values = None

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
