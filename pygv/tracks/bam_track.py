import os
from collections import namedtuple

import numpy as np
import pysam
from matplotlib.collections import PatchCollection
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from pygv.errors.DataIntegrity import BamIndexDoesntExists

from .bed_track import _LaneRegistry
from .track import NumericalTrack, Track


class _GenericNumericalBamTrack(NumericalTrack):
    """
    Generic numerical track for bam files

    Parameters
    ----------
    track : str
        Path to the bam file
    filters : None or a callable function
        :attr:`filters`
    kwargs :
        read_colors : list of color-like values
            :attr:`read_colors`
        flip_strand : bool
            :attr:`flip_strand`
    """

    @property
    def read_colors(self):
        """
        Colors for reads in different conditions, like forward/reverse. Default `("#E69696", "#9696E6")`
        """
        return self._read_colors

    @read_colors.setter
    def read_colors(self, value):
        if all(map(is_color_like, value)):
            self._read_colors = value
        else:
            self._read_colors = ("#E69696", "#9696E6")

    @property
    def flip_strand(self):
        """
        Flip the strand of reads, by default `False`
        """
        return self._flip_strand

    @flip_strand.setter
    def flip_strand(self, value):
        try:
            self._flip_strand = bool(value)
        except:
            pass

    @property
    def filters(self):
        """
        None or a function which returns True/False for each read, reads with Trues will be kept
        By default `None`
        """
        return self._filters

    @filters.setter
    def filters(self, value):
        if value is None or callable(value):
            self._filters = value
        else:
            print("Invalid filter")

    def __init__(self, track, filters=None, **kwargs):
        super(_GenericNumericalBamTrack, self).__init__(**kwargs)
        if not os.path.exists(track):
            raise ValueError

        self.bam = pysam.AlignmentFile(track)
        if not self.bam.has_index():
            raise BamIndexDoesntExists(
                "Cannot locate index for {bam}".format(bam=track)
            )

        self._read_colors = ("#E69696", "#9696E6")
        self.read_colors = kwargs.pop("read_colors", ("#E69696", "#9696E6"))
        self._flip_strand = False
        self.flip_strand = kwargs.pop("flip_strand", False)
        self._filters = None
        self.filters = filters

    def _get(self, chromosome, start, end):
        pass


class _GenericBamTrack(Track):
    """
    Generic track for bam files

    Parameters
    ----------
    track : str
        Path to the bam file
    filters : None or a callable function
        :attr:`filters`
    kwargs :
        color_reads_by : str or None
            :attr:`color_reads_by`
        color_legends : None or list
            :attr:`color_legends`
        read_colors : list of color-like values
            :attr:`read_colors`
        flip_strand : bool
            :attr:`flip_strand`
    """

    @property
    def read_colors(self):
        """
        Colors for reads in different conditions, like forward/reverse. Default `("#E69696", "#9696E6")`
        """
        return self._read_colors

    @read_colors.setter
    def read_colors(self, value):
        if all(map(is_color_like, value)):
            self._read_colors = value
        else:
            self._read_colors = ("#E69696", "#9696E6")

    @property
    def color_reads_by(self):
        """
        Color reads by certain criteria, currently supported values:
        * first of pair strand
        * read strand
        set it to None to disable this function
        """
        return self._color_reads_by

    @color_reads_by.setter
    def color_reads_by(self, value):
        if value in self._SUPPORTED_COLOR_METHOD or value is None:
            self._color_reads_by = value
        else:
            raise ValueError("Not supported option")

    @property
    def filters(self):
        """
        None or a function which returns True/False for each read, reads with Trues will be kept
        By default `None`
        """
        return self._filters

    @filters.setter
    def filters(self, value):
        if value is None or callable(value):
            self._filters = value
        else:
            print("Invalid filter")

    @property
    def color_legends(self):
        """
        List of strs for each legend, set it as None to disable legend
        """
        return self._color_legends

    @color_legends.setter
    def color_legends(self, value):
        try:
            self._color_legends = value
        except:
            pass

    @property
    def legend_title(self):
        """
        List of strs for each legend, set it as None to disable legend
        """
        return self._legend_title

    @legend_title.setter
    def legend_title(self, value):
        try:
            self._legend_title = value
        except:
            pass

    @property
    def sampling_ratio(self):
        return self._sampling_ratio

    @sampling_ratio.setter
    def sampling_ratio(self, value):
        self._sampling_ratio = float(value)

    def __init__(self, track, filters=None, **kwargs):
        super(_GenericBamTrack, self).__init__(**kwargs)
        if not os.path.exists(track):
            raise ValueError

        self.bam = pysam.AlignmentFile(track)
        if not self.bam.has_index():
            raise BamIndexDoesntExists(
                "Cannot locate index for {bam}".format(bam=track)
            )

        # lane manager
        self._lane_registries = []
        self.allowed_features = kwargs.pop("allowed_features", None)

        # read colors
        self._SUPPORTED_COLOR_METHOD = (
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
        )
        self._color_reads_by = None
        self.color_reads_by = kwargs.pop("color_reads_by", None)
        self._color_legends = None
        self.color_legends = kwargs.pop("color_legends", None)
        self._legend_title = None
        self.legend_title = kwargs.pop("legend_title", "Mapping direction")

        self._read_colors = ("#E69696", "#9696E6")
        self.read_colors = kwargs.pop("read_colors", ("#E69696", "#9696E6"))

        self._allowed_features = None
        self.allowed_features = kwargs.pop("allowed_features", None)

        self._filters = None
        self.filters = filters

        self._sampling_ratio = 1.0

    def _get(self, chromosome, start, end):
        pass

    def _assign_read_color(self, read):
        color = self._read_colors[0]
        if self._color_reads_by == "first of pair strand":
            if read.is_read1:
                if read.is_reverse:
                    color = self._read_colors[1]
        elif self._color_reads_by == "read strand":
            if read.is_reverse:
                color = self._read_colors[1]
        return color


class CoverageTrack(_GenericNumericalBamTrack):
    """
    Bam coverage track

    Examples
    --------

    .. plot:: ../examples/plot_bam_coverage.py
    """

    def _get(self, chromosome, start, end):
        func = "all" if self._filters is None else self._filters
        a, c, g, t = self.bam.count_coverage(chromosome, start, end, read_callback=func)
        cov = []
        for i, j, k, l in zip(a, c, g, t):
            cov.append(i + j + k + l)
        values = np.array(cov)
        xvalues = np.arange(start, end, step=1)

        if self._stat_method is not None:
            from scipy.stats import binned_statistic

            y_new, x_new, _ = binned_statistic(
                xvalues, values, statistic=self._stat_method, bins=self._n_bins
            )
            xvalues = x_new
            values = y_new
        keep_idx = self._merge_redundant_values(xvalues, values)
        return xvalues[keep_idx], values[keep_idx]

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw coverage track from bam file

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the ROI/window, 0-based
        end : int
            end of the ROI/window, 0-based
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(NumericalTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        x, y = self._get(chromosome=chromosome, start=start, end=end)
        self._ax.plot(x, y, color=self.color, linewidth=self.line_width)
        # self.ax.bar(x, y, color=self.color, width=1)
        self._ax.fill_between(
            x, y, 0, facecolor=self.color, alpha=self._alpha, lw=self.line_width
        )
        try:
            self._ax.ticklabel_format(style="plain", useOffset=False)
        except:
            pass


class CollapsedReadTrack(_GenericBamTrack):
    """
    Plot collapsed reads (only 5' end and the span)

    Parameters
    ----------
    track :
    kwargs :
        patch_height :
            :attr:`patch_height`
        line_color :
            :attr:`line_color`
        max_num_read :
            :attr:`max_num_read`
        pileup_offset :
            :attr:`pileup_offset`

    Examples
    --------

    .. plot:: ../examples/plot_bam_collapsed_reads.py
    """

    @property
    def line_color(self):
        """
        Line color, by default, #96B8C8
        """
        return self._line_color

    @line_color.setter
    def line_color(self, value):
        try:
            self._line_color = value
        except:
            pass

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
    def max_num_read(self):
        """
        Max number of reads in the visible window, if there are more reads, random downsampling will be used
        """
        return self._patch_height

    @max_num_read.setter
    def max_num_read(self, value):
        try:
            self._max_num_read = int(value)
        except:
            pass

    @property
    def pileup_offset(self):
        """
        offset of pileup, by default, 0.1
        """
        return self.pileup_offset

    @pileup_offset.setter
    def pileup_offset(self, value):
        try:
            self._pileup_offset = float(value)
        except:
            pass

    def __init__(self, track, **kwargs):
        super(CollapsedReadTrack, self).__init__(track, **kwargs)
        self._patch_height = 1
        self.patch_height = kwargs.pop("patch_height", 1)
        self._line_color = "#96B8C8"
        self.line_color = kwargs.pop("line_color", "#96B8C8")
        self._max_num_read = 500
        self.max_num_read = kwargs.pop("max_num_read", 500)
        self._pileup_offset = 0.1
        self.pileup_offset = kwargs.pop("pileup_offset", 0.1)

    def _get(self, chromosome, start, end):
        # values = np.nan_to_num(self.bw.values(chromosome, start, end))
        func = "all" if self._filters is None else self._filters
        reads = []
        Read = namedtuple(
            "Read", field_names=("start", "end", "length", "strand", "aligned_segment")
        )
        for read in self.bam.fetch(contig=chromosome, start=start, stop=end):
            if func != "all":
                if not func(read):
                    continue
            if np.random.random() > self.sampling_ratio:
                continue
            reads.append(
                Read._make(
                    (
                        read.reference_start,
                        read.reference_end,
                        read.reference_length,
                        -1 if read.is_reverse else 1,
                        read,
                    )
                )
            )
        return reads

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        pass

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw collapsed-read track

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the ROI/window, 0-based
        end : int
            end of the ROI/window, 0-based
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(CollapsedReadTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )

        read_records = []
        read_records_dict = {}
        for read in self._get(chromosome=chromosome, start=start, end=end):
            k = "%d-%d" % (read.start, read.length)
            read_records.append((read.start, read.length))
            if k not in read_records_dict.keys():
                read_records_dict[k] = [read.start, read.length, 1]
            else:
                read_records_dict[k][2] += 1

        rrs = sorted(read_records, key=lambda x: x[0])
        Y = np.linspace(0, 1, len(read_records))
        sr = 1
        if len(rrs) > self._max_num_read:
            sr = self._max_num_read / len(rrs)
        for i, rr in enumerate(rrs):
            if np.random.random() < sr:
                y_offset = (
                    np.random.random() * self._pileup_offset - self._pileup_offset
                )
                self._ax.plot(
                    (rr[0], rr[0] + rr[1]),
                    (Y[i] + y_offset, Y[i] + y_offset),
                    zorder=1,
                    alpha=self._alpha / 2,
                    color=self._line_color,
                    lw=self.line_width,
                )
                self._ax.scatter(
                    rr[0],
                    Y[i] + y_offset,
                    marker="|",
                    color=self.color,
                    alpha=self._alpha,
                    zorder=2,
                    s=5,
                )

        self._ax.yaxis.set_ticks([])
        self._ax.set_xlim((start, end))


class SplicedReadTrack(_GenericBamTrack):
    """
    Plot spliced reads

    Parameters
    ----------
    track :
    kwargs :
        padding_left :
            :attr:`padding_left`
        padding_right :
            :attr:`padding_right`
        show_name :
            :attr:`show_name`
        patch_height :
            :attr:`patch_height`
        lane_space :
            :attr:`lane_space`
        features_per_lane :
            :attr:`features_per_lane`
        line_color :
            :attr:`line_color`

    Examples
    --------

    .. plot:: ../examples/plot_bam_spliced_reads.py
    """

    @property
    def padding_left(self):
        """
        Units adding to the left of features (adding places for text labels)
        """
        return self._padding_left

    @padding_left.setter
    def padding_left(self, value):
        try:
            self._padding_left = float(value)
        except:
            pass

    @property
    def padding_right(self):
        """
        Units adding to the right of features (adding places for text labels)
        """
        return self._padding_right

    @padding_right.setter
    def padding_right(self, value):
        try:
            self._padding_right = float(value)
        except:
            pass

    @property
    def show_name(self):
        """
        Units adding to the right of features (adding places for text labels)
        """
        return self._show_name

    @show_name.setter
    def show_name(self, value):
        try:
            self._show_name = bool(value)
        except:
            pass

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

    def __init__(self, track, **kwargs):
        super(SplicedReadTrack, self).__init__(track, **kwargs)

        self._lane_registries = []
        self._padding_left = 0
        self.padding_left = kwargs.pop("padding_left", 0)
        self._padding_right = 0
        self.padding_right = kwargs.pop("padding_right", 0)
        self._show_name = True
        self.show_name = kwargs.pop("show_name", True)
        self._patch_height = 1
        self.patch_height = kwargs.pop("patch_height", 1)
        self._lane_space = 0.25
        self.lane_space = kwargs.pop("lane_space", 0.25)
        self._line_color = "black"
        self.line_color = kwargs.pop("line_color", "black")
        self._box_color = kwargs.pop("box_color", "#A1A1A1")
        self._box_border = kwargs.pop("box_border", "#6E6E6E")
        self._features_per_lane = 3
        self.features_per_lane = kwargs.pop("features_per_lane", 3)

    def _get(self, chromosome, start, end):
        func = "all" if self._filters is None else self._filters
        reads = []
        Read = namedtuple(
            "Read", field_names=("start", "end", "length", "strand", "aligned_segment")
        )
        for read in self.bam.fetch(contig=chromosome, start=start, stop=end):
            if func != "all":
                if not func(read):
                    continue
            if np.random.random() > self.sampling_ratio:
                continue
            reads.append(
                Read._make(
                    (
                        read.reference_start,
                        read.reference_end,
                        read.reference_length,
                        -1 if read.is_reverse else 1,
                        read,
                    )
                )
            )
        return reads

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        """
        Build non-overlapping tracks

        Parameters
        ----------
        chromosome : str
            Chromosome
        start : int
            start of visible window
        end : int
            end of visible window
        Returns
        -------

        """
        # clean up lane registry
        self._lane_registries = []
        for interval in self._get(chromosome=chromosome, start=start, end=end):
            if self.color_reads_by == "first of pair strand":
                if not interval.aligned_segment.is_read1:
                    continue

            active_lane = None
            start_loc = int(interval.start)
            end_loc = int(interval.end)

            if len(self._lane_registries) == 0:
                self._lane_registries.append(_LaneRegistry())

            for lr in self._lane_registries:
                if lr.max_coord is None:
                    active_lane = lr.offset
                    lr.min_coord = start_loc
                    lr.max_coord = end_loc
                    lr.features.append(interval)
                    break
                else:
                    if lr.max_coord < start_loc - self._padding_left:
                        active_lane = lr.offset
                        lr.min_coord = min(lr.min_coord, start_loc)
                        lr.max_coord = max(lr.max_coord, end_loc)
                        lr.features.append(interval)
                        break

            if (
                type(self._allowed_features) is int
                and len(self._lane_registries) >= self._allowed_features
                and active_lane is None
            ):
                continue

            if active_lane is None:
                self._lane_registries.append(
                    _LaneRegistry(
                        offset=len(self._lane_registries),
                        min_coord=start_loc,
                        max_coord=end_loc,
                        features=[interval],
                    )
                )

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw spliced-read track

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the ROI/window, 0-based
        end : int
            end of the ROI/window, 0-based
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(SplicedReadTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        for lane in self._lane_registries:
            for interval in lane.features:
                start_loc = int(interval.start)
                end_loc = int(interval.end)
                visible_start = max(start_loc, start)
                visible_end = min(end_loc, end)
                active_lane = lane.offset

                real_active_line = (self._patch_height + self._lane_space) * active_lane

                blocks = interval.aligned_segment.get_blocks()
                self._ax.plot(
                    (
                        start_loc if start_loc >= start else visible_start,
                        end_loc - 1 if end_loc <= end else visible_end,
                    ),
                    (-1 * real_active_line, -1 * real_active_line),
                    # color=self._line_color,
                    color=self._assign_read_color(interval.aligned_segment),
                    alpha=self._alpha,
                    linewidth=self.line_width,
                    clip_on=True,
                )

                patches = []
                patch_colors = []
                for block in blocks:
                    x = block[0]
                    if x < end or block[1] > start:
                        # in case overflow
                        # adjusted_size = size
                        plot_start = (
                            block[0] if block[0] >= visible_start else visible_start
                        )
                        plot_end = block[1] if block[1] <= visible_end else visible_end
                        p = Rectangle(
                            xy=(
                                plot_start,
                                -1 * real_active_line - (self._patch_height / 2),
                            ),
                            width=plot_end - plot_start,
                            clip_on=True,
                            height=self._patch_height,
                        )
                        patches.append(p)
                        patch_colors.append(
                            self._assign_read_color(interval.aligned_segment)
                        )

                self._ax.add_collection(
                    PatchCollection(
                        patches,
                        edgecolors=patch_colors,
                        facecolors=patch_colors,
                        linewidths=self.line_width,
                        zorder=20,
                        clip_on=True,
                    )
                )

        if self._color_legends is not None:
            legend_lines = []
            if len(self._color_legends) == len(self._read_colors):
                for c in self._read_colors:
                    legend_lines.append(Line2D([0], [0], color=c, lw=1))
            self._ax.legend(
                legend_lines,
                self._color_legends,
                title=self._legend_title,
                frameon=False,
            )
        self._ax.yaxis.set_ticks([])
        self._ax.set_xlim((start, end))


class StrandSpecificCoverageTrack(_GenericNumericalBamTrack):
    """
    Draw strand-specific coverages from bam file

    Examples
    --------

    .. plot:: ../examples/plot_bam_stranded_coverage.py
    """

    def _get(self, chromosome, start, end, direction="forward"):
        if direction == "forward":
            if self._filters is None:
                if not self._flip_strand:
                    func = lambda x: not x.is_reverse
                else:
                    func = lambda x: x.is_reverse
            else:
                if not self._flip_strand:
                    func = lambda x: not x.is_reverse and self._filters(x)
                else:
                    func = lambda x: x.is_reverse and self._filters(x)
        else:
            if self._filters is None:
                if not self._flip_strand:
                    func = lambda x: x.is_reverse
                else:
                    func = lambda x: not x.is_reverse
            else:
                if not self._flip_strand:
                    func = lambda x: x.is_reverse and self._filters(x)
                else:
                    func = lambda x: not x.is_reverse and self._filters(x)

        a, c, g, t = self.bam.count_coverage(chromosome, start, end, read_callback=func)
        cov = []
        for i, j, k, l in zip(a, c, g, t):
            cov.append(i + j + k + l)
        values = np.array(cov)

        if self.scale != 1:
            values = self.scale * values

        xvalues = np.arange(start, end, step=1)

        if self.stat_method is not None:
            from scipy.stats import binned_statistic

            y_new, x_new, _ = binned_statistic(
                xvalues, values, statistic=self.stat_method, bins=self.n_bins
            )
            xvalues = x_new
            values = y_new
        keep_idx = self._merge_redundant_values(xvalues, values)
        return xvalues[keep_idx], values[keep_idx]

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw strand-specific coverages from bam file

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the ROI/window, 0-based
        end : int
            end of the ROI/window, 0-based
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(NumericalTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        x, yp = self._get(chromosome=chromosome, start=start, end=end)
        self._ax.plot(x, yp, color=self._read_colors[0], linewidth=self.line_width)
        self._ax.fill_between(
            x, yp, 0, facecolor=self._read_colors[0], alpha=self._alpha
        )
        x, ym = self._get(
            chromosome=chromosome, start=start, end=end, direction="reverse"
        )
        self._ax.plot(x, -1 * ym, color=self._read_colors[1], linewidth=self.line_width)
        self._ax.fill_between(
            x, -1 * ym, 0, facecolor=self._read_colors[1], alpha=self._alpha
        )
        # self.ax.bar(x, y, color=self.color, width=1)

        try:
            self._ax.ticklabel_format(style="plain", useOffset=False)
        except:
            pass


class ReadArcTrack(_GenericBamTrack):
    """
    Plot read in arcs

    Examples
    --------

    .. plot:: ../examples/plot_bam_arc_reads.py
    """

    def _get(self, chromosome, start, end):
        func = "all" if self._filters is None else self._filters
        reads = []
        Read = namedtuple("Read", field_names=("start", "end", "length", "strand"))
        for read in self.bam.fetch(contig=chromosome, start=start, stop=end):
            if np.random.random() > self.sampling_ratio:
                continue
            reads.append(
                Read._make(
                    (
                        read.reference_start,
                        read.reference_end,
                        read.reference_length,
                        -1 if read.is_reverse else 1,
                    )
                )
            )
        return reads

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super(ReadArcTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        import matplotlib as mpl
        from matplotlib.patches import Arc

        reads = self._get(chromosome=chromosome, start=start, end=end)
        starts = []
        ends = []
        for read in reads:
            x = read.start + read.length / 2 * read.strand
            # linewidth=cnt / total * linewidth_factor
            arc = Arc(
                (x, 0),
                read.length,
                read.length / 2 * read.strand,
                theta2=180,
                color=self.color,
                alpha=self.alpha,
            )
            ax.add_patch(arc)
            starts.append(read.start)
            ends.append(read.end)

        start_color = kwargs.get("start_color", "green")
        end_color = kwargs.get("start_color", "red")
        self._ax.scatter(
            starts,
            [
                0,
            ]
            * len(starts),
            color=start_color,
            alpha=self._alpha,
            s=self.line_width,
        )
        self._ax.scatter(
            ends,
            [
                0,
            ]
            * len(ends),
            color=end_color,
            alpha=self._alpha,
            s=self.line_width,
        )
        self._ax.yaxis.set_ticks([])
        self._ax.set_xlim((start, end))
        # self._ax.ticklabel_format(style="plain", useOffset=False)
        self._ax.autoscale_view()
