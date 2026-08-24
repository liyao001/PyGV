import os
from collections import namedtuple
from typing import Any, Optional, Union

import numpy as np
import pysam
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pydantic import Field, PrivateAttr

from pygv.errors.DataIntegrity import BamIndexDoesntExists

from .bed_track import _LaneRegistry
from .track import NumericalTrack, Track
from .types import Color, ColorReadsBy, ColorSequence, FilterFn


class _GenericNumericalBamTrack(NumericalTrack):
    """Generic numerical track for BAM files."""

    track: str = Field(description="Path to the BAM file")
    filters: FilterFn = Field(
        default=None,
        description="Callable returning True/False for each read; None keeps all reads",
    )
    read_colors: ColorSequence = Field(
        default=("#E69696", "#9696E6"),
        kw_only=True,
        description="Colors for reads in different conditions, like forward/reverse",
    )
    flip_strand: bool = Field(
        default=False, kw_only=True, description="Flip the strand of reads"
    )

    _bam: Any = PrivateAttr(default=None)

    def __init__(self, track: str, filters: FilterFn = None, **data: Any) -> None:
        super().__init__(track=track, filters=filters, **data)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        if not os.path.exists(self.track):
            raise ValueError
        self._bam = pysam.AlignmentFile(self.track)
        if not self._bam.has_index():
            raise BamIndexDoesntExists(
                "Cannot locate index for {bam}".format(bam=self.track)
            )

    def _get(self, chromosome, start, end):
        pass


class _GenericBamTrack(Track):
    """Generic track for BAM files."""

    track: str = Field(description="Path to the BAM file")
    filters: FilterFn = Field(
        default=None,
        description="Callable returning True/False for each read; None keeps all reads",
    )
    allowed_features: Optional[int] = Field(
        default=None,
        kw_only=True,
        description="Max number of feature lanes to plot",
    )
    color_reads_by: Optional[ColorReadsBy] = Field(
        default=None,
        kw_only=True,
        description="Color reads by a supported criterion, or None to disable",
    )
    color_legends: Optional[Union[list, tuple]] = Field(
        default=None,
        kw_only=True,
        description="Legend labels; None disables the legend",
    )
    legend_title: Optional[str] = Field(
        default="Mapping direction",
        kw_only=True,
        description="Title for the read-color legend",
    )
    read_colors: ColorSequence = Field(
        default=("#E69696", "#9696E6"),
        kw_only=True,
        description="Colors for reads in different conditions, like forward/reverse",
    )
    sampling_ratio: float = Field(
        default=1.0, kw_only=True, description="Fraction of reads to sample"
    )

    _bam: Any = PrivateAttr(default=None)
    _lane_registries: list = PrivateAttr(default_factory=list)

    def __init__(self, track: str, filters: FilterFn = None, **data: Any) -> None:
        super().__init__(track=track, filters=filters, **data)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        if not os.path.exists(self.track):
            raise ValueError
        self._bam = pysam.AlignmentFile(self.track)
        if not self._bam.has_index():
            raise BamIndexDoesntExists(
                "Cannot locate index for {bam}".format(bam=self.track)
            )

    def _get(self, chromosome, start, end):
        pass

    def _assign_read_color(self, read):
        color = self.read_colors[0]
        if self.color_reads_by == "first of pair strand":
            if read.is_read1:
                if read.is_reverse:
                    color = self.read_colors[1]
        elif self.color_reads_by == "read strand":
            if read.is_reverse:
                color = self.read_colors[1]
        return color


class CoverageTrack(_GenericNumericalBamTrack):
    """
    Bam coverage track

    Examples
    --------

    .. plot:: ../examples/plot_bam_coverage.py
    """

    def _get(self, chromosome, start, end):
        func = "all" if self.filters is None else self.filters
        a, c, g, t = self._bam.count_coverage(chromosome, start, end, read_callback=func)
        cov = []
        for i, j, k, l in zip(a, c, g, t):
            cov.append(i + j + k + l)
        values = np.array(cov)
        values = self.data_transform(values)
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
            x, y, 0, facecolor=self.color, alpha=self.alpha, lw=self.line_width
        )
        try:
            self._ax.ticklabel_format(style="plain", useOffset=False)
        except:
            pass


class CollapsedReadTrack(_GenericBamTrack):
    """Plot collapsed reads (only 5' end and the span)."""

    patch_height: float = Field(default=1, kw_only=True)
    line_color: Color = Field(default="#96B8C8", kw_only=True)
    max_num_read: int = Field(
        default=500,
        kw_only=True,
        description="Max number of reads in the window; extra reads are downsampled",
    )
    pileup_offset: float = Field(default=0.1, kw_only=True)

    def _get(self, chromosome, start, end):
        # values = np.nan_to_num(self.bw.values(chromosome, start, end))
        func = "all" if self.filters is None else self.filters
        reads = []
        Read = namedtuple(
            "Read", field_names=("start", "end", "length", "strand", "aligned_segment")
        )
        for read in self._bam.fetch(contig=chromosome, start=start, stop=end):
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
        if len(rrs) > self.max_num_read:
            sr = self.max_num_read / len(rrs)
        for i, rr in enumerate(rrs):
            if np.random.random() < sr:
                y_offset = (
                    np.random.random() * self.pileup_offset - self.pileup_offset
                )
                self._ax.plot(
                    (rr[0], rr[0] + rr[1]),
                    (Y[i] + y_offset, Y[i] + y_offset),
                    zorder=1,
                    alpha=self.alpha / 2,
                    color=self.line_color,
                    lw=self.line_width,
                )
                self._ax.scatter(
                    rr[0],
                    Y[i] + y_offset,
                    marker="|",
                    color=self.color,
                    alpha=self.alpha,
                    zorder=2,
                    s=5,
                )

        self._ax.yaxis.set_ticks([])
        self._ax.set_xlim((start, end))


class SplicedReadTrack(_GenericBamTrack):
    """Plot spliced reads."""

    padding_left: float = Field(default=0, kw_only=True)
    padding_right: float = Field(default=0, kw_only=True)
    show_name: bool = Field(default=True, kw_only=True)
    patch_height: float = Field(default=1, kw_only=True)
    lane_space: float = Field(default=0.25, kw_only=True)
    line_color: Color = Field(default="black", kw_only=True)
    box_color: Color = Field(default="#A1A1A1", kw_only=True)
    box_border: Color = Field(default="#6E6E6E", kw_only=True)
    features_per_lane: int = Field(default=3, kw_only=True)

    def _get(self, chromosome, start, end):
        func = "all" if self.filters is None else self.filters
        reads = []
        Read = namedtuple(
            "Read", field_names=("start", "end", "length", "strand", "aligned_segment")
        )
        for read in self._bam.fetch(contig=chromosome, start=start, stop=end):
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
                    if lr.max_coord < start_loc - self.padding_left:
                        active_lane = lr.offset
                        lr.min_coord = min(lr.min_coord, start_loc)
                        lr.max_coord = max(lr.max_coord, end_loc)
                        lr.features.append(interval)
                        break

            if (
                type(self.allowed_features) is int
                and len(self._lane_registries) >= self.allowed_features
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

                real_active_line = (self.patch_height + self.lane_space) * active_lane

                blocks = interval.aligned_segment.get_blocks()
                self._ax.plot(
                    (
                        start_loc if start_loc >= start else visible_start,
                        end_loc - 1 if end_loc <= end else visible_end,
                    ),
                    (-1 * real_active_line, -1 * real_active_line),
                    # color=self.line_color,
                    color=self._assign_read_color(interval.aligned_segment),
                    alpha=self.alpha,
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
                                -1 * real_active_line - (self.patch_height / 2),
                            ),
                            width=plot_end - plot_start,
                            clip_on=True,
                            height=self.patch_height,
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

        if self.color_legends is not None:
            legend_lines = []
            if len(self.color_legends) == len(self.read_colors):
                for c in self.read_colors:
                    legend_lines.append(Line2D([0], [0], color=c, lw=1))
            self._ax.legend(
                legend_lines,
                self.color_legends,
                title=self.legend_title,
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
            if self.filters is None:
                if not self.flip_strand:
                    func = lambda x: not x.is_reverse
                else:
                    func = lambda x: x.is_reverse
            else:
                if not self.flip_strand:
                    func = lambda x: not x.is_reverse and self.filters(x)
                else:
                    func = lambda x: x.is_reverse and self.filters(x)
        else:
            if self.filters is None:
                if not self.flip_strand:
                    func = lambda x: x.is_reverse
                else:
                    func = lambda x: not x.is_reverse
            else:
                if not self.flip_strand:
                    func = lambda x: x.is_reverse and self.filters(x)
                else:
                    func = lambda x: not x.is_reverse and self.filters(x)

        a, c, g, t = self._bam.count_coverage(chromosome, start, end, read_callback=func)
        cov = []
        for i, j, k, l in zip(a, c, g, t):
            cov.append(i + j + k + l)
        values = np.array(cov)
        values = self.data_transform(values)
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
        self._ax.plot(x, yp, color=self.read_colors[0], linewidth=self.line_width)
        self._ax.fill_between(
            x, yp, 0, facecolor=self.read_colors[0], alpha=self.alpha
        )
        x, ym = self._get(
            chromosome=chromosome, start=start, end=end, direction="reverse"
        )
        self._ax.plot(x, -1 * ym, color=self.read_colors[1], linewidth=self.line_width)
        self._ax.fill_between(
            x, -1 * ym, 0, facecolor=self.read_colors[1], alpha=self.alpha
        )
        # self.ax.bar(x, y, color=self.color, width=1)

        try:
            self._ax.ticklabel_format(style="plain", useOffset=False)
        except:
            pass


class ReadArcTrack(_GenericBamTrack):
    """Plot reads as arcs."""

    start_color: Color = Field(default="green", kw_only=True)
    end_color: Color = Field(default="red", kw_only=True)

    def _get(self, chromosome, start, end):
        func = "all" if self.filters is None else self.filters
        reads = []
        Read = namedtuple("Read", field_names=("start", "end", "length", "strand"))
        for read in self._bam.fetch(contig=chromosome, start=start, stop=end):
            if func != "all" and not func(read):
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
                    )
                )
            )
        return reads

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super(ReadArcTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
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

        start_color = self.start_color
        end_color = self.end_color
        self._ax.scatter(
            starts,
            [
                0,
            ]
            * len(starts),
            color=start_color,
            alpha=self.alpha,
            s=self.line_width,
        )
        self._ax.scatter(
            ends,
            [
                0,
            ]
            * len(ends),
            color=end_color,
            alpha=self.alpha,
            s=self.line_width,
        )
        self._ax.yaxis.set_ticks([])
        self._ax.set_xlim((start, end))
        # self._ax.ticklabel_format(style="plain", useOffset=False)
        self._ax.autoscale_view()
