import os
from collections import namedtuple
from typing import ClassVar

import numpy as np
import pysam
from matplotlib.collections import PatchCollection
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pydantic import Field, field_validator

from pygv.errors.DataIntegrity import BamIndexDoesntExists

from .bed_track import _LaneRegistry
from .track import NumericalTrack, NumericalTrackConfig, Track, TrackConfig


class GenericNumericalBamTrackConfig(NumericalTrackConfig):
    read_colors: object = Field(
        default=("#E69696", "#9696E6"),
        description="Matplotlib colors used for reads in different strand conditions.",
    )
    flip_strand: bool = Field(
        default=False, description="Whether to flip read strand orientation."
    )

    @field_validator("read_colors")
    def _validate_read_colors(cls, value):
        if all(map(is_color_like, value)):
            return value
        raise ValueError("read_colors must contain only Matplotlib color-like values")


class GenericBamTrackConfig(TrackConfig):
    allowed_features: object = None
    color_reads_by: object = Field(
        default=None, description="Read attribute or method used to color alignments."
    )
    color_legends: object = Field(
        default=None, description="Legend labels for read color categories."
    )
    legend_title: object = Field(
        default="Mapping direction", description="Title shown above the read color legend."
    )
    read_colors: object = Field(
        default=("#E69696", "#9696E6"),
        description="Matplotlib colors used for read categories.",
    )
    sampling_ratio: float = Field(
        default=1.0, ge=0, le=1, description="Fraction of reads to sample for plotting."
    )

    _SUPPORTED_COLOR_METHOD: ClassVar[tuple] = (
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

    @field_validator("read_colors")
    def _validate_read_colors(cls, value):
        if all(map(is_color_like, value)):
            return value
        raise ValueError("read_colors must contain only Matplotlib color-like values")

    @field_validator("color_reads_by")
    def _validate_color_reads_by(cls, value):
        if value in cls._SUPPORTED_COLOR_METHOD or value is None:
            return value
        raise ValueError("Not supported option")


class CollapsedReadTrackConfig(GenericBamTrackConfig):
    patch_height: float = Field(
        default=1, gt=0, description="Height of each collapsed read patch."
    )
    line_color: object = Field(
        default="#96B8C8", description="Matplotlib color used for collapsed read lines."
    )
    max_num_read: int = Field(
        default=500, gt=0, description="Maximum number of reads to draw."
    )
    pileup_offset: float = Field(
        default=0.1, ge=0, description="Vertical offset between piled-up reads."
    )

    @field_validator("line_color")
    def _validate_line_color(cls, value):
        if value is None or is_color_like(value):
            return value
        raise ValueError(f"Invalid color value: {value}")


class SplicedReadTrackConfig(GenericBamTrackConfig):
    padding_left: float = Field(
        default=0, ge=0, description="Extra left-side spacing used when placing reads."
    )
    padding_right: float = Field(
        default=0, ge=0, description="Extra right-side spacing used when placing reads."
    )
    show_name: bool = Field(default=True, description="Whether to draw read names.")
    patch_height: float = Field(
        default=1, gt=0, description="Height of each spliced read patch."
    )
    lane_space: float = Field(
        default=0.25, ge=0, description="Extra vertical spacing between read lanes."
    )
    features_per_lane: int = Field(
        default=3, gt=0, description="Maximum number of reads grouped into each lane."
    )
    line_color: object = Field(
        default="black", description="Matplotlib color used for spliced read connector lines."
    )
    box_color: object = Field(
        default="#A1A1A1", description="Matplotlib fill color for spliced read boxes."
    )
    box_border: object = Field(
        default="#6E6E6E", description="Matplotlib edge color for spliced read boxes."
    )

    @field_validator("line_color", "box_color", "box_border")
    def _validate_colors(cls, value):
        if value is None or is_color_like(value):
            return value
        raise ValueError(f"Invalid color value: {value}")


class _GenericNumericalBamTrack(NumericalTrack, GenericNumericalBamTrackConfig):
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

    _FIELD_PRIVATE_ATTRS: ClassVar[dict] = {
        **NumericalTrack._FIELD_PRIVATE_ATTRS,
        "read_colors": "_read_colors",
        "flip_strand": "_flip_strand",
    }

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
        self.read_colors = self.__dict__["read_colors"]
        self._flip_strand = False
        self.flip_strand = self.__dict__["flip_strand"]
        self._filters = None
        self.filters = filters

    def _get(self, chromosome, start, end):
        pass


class _GenericBamTrack(Track, GenericBamTrackConfig):
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

    _FIELD_PRIVATE_ATTRS: ClassVar[dict] = {
        **Track._FIELD_PRIVATE_ATTRS,
        "allowed_features": "_allowed_features",
        "color_reads_by": "_color_reads_by",
        "color_legends": "_color_legends",
        "legend_title": "_legend_title",
        "read_colors": "_read_colors",
        "sampling_ratio": "_sampling_ratio",
    }

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
        self.allowed_features = self.__dict__["allowed_features"]

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
        self.color_reads_by = self.__dict__["color_reads_by"]
        self._color_legends = None
        self.color_legends = self.__dict__["color_legends"]
        self._legend_title = None
        self.legend_title = self.__dict__["legend_title"]

        self._read_colors = ("#E69696", "#9696E6")
        self.read_colors = self.__dict__["read_colors"]

        self._allowed_features = None
        self.allowed_features = self.__dict__["allowed_features"]

        self._filters = None
        self.filters = filters

        self._sampling_ratio = self.__dict__["sampling_ratio"]

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


class CollapsedReadTrack(_GenericBamTrack, CollapsedReadTrackConfig):
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

    _FIELD_PRIVATE_ATTRS: ClassVar[dict] = {
        **_GenericBamTrack._FIELD_PRIVATE_ATTRS,
        "patch_height": "_patch_height",
        "line_color": "_line_color",
        "max_num_read": "_max_num_read",
        "pileup_offset": "_pileup_offset",
    }

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
        return self._max_num_read

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
        return self._pileup_offset

    @pileup_offset.setter
    def pileup_offset(self, value):
        try:
            self._pileup_offset = float(value)
        except:
            pass

    def __init__(self, track, **kwargs):
        config = CollapsedReadTrackConfig(**kwargs)
        super(CollapsedReadTrack, self).__init__(track, **kwargs)
        self._patch_height = 1
        self.patch_height = config.patch_height
        self._line_color = "#96B8C8"
        self.line_color = config.line_color
        self._max_num_read = 500
        self.max_num_read = config.max_num_read
        self._pileup_offset = 0.1
        self.pileup_offset = config.pileup_offset

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


class SplicedReadTrack(_GenericBamTrack, SplicedReadTrackConfig):
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

    _FIELD_PRIVATE_ATTRS: ClassVar[dict] = {
        **_GenericBamTrack._FIELD_PRIVATE_ATTRS,
        "padding_left": "_padding_left",
        "padding_right": "_padding_right",
        "show_name": "_show_name",
        "patch_height": "_patch_height",
        "lane_space": "_lane_space",
        "features_per_lane": "_features_per_lane",
        "line_color": "_line_color",
        "box_color": "_box_color",
        "box_border": "_box_border",
    }

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
        config = SplicedReadTrackConfig(**kwargs)
        super(SplicedReadTrack, self).__init__(track, **kwargs)

        self._lane_registries = []
        self._padding_left = 0
        self.padding_left = config.padding_left
        self._padding_right = 0
        self.padding_right = config.padding_right
        self._show_name = True
        self.show_name = config.show_name
        self._patch_height = 1
        self.patch_height = config.patch_height
        self._lane_space = 0.25
        self.lane_space = config.lane_space
        self._line_color = "black"
        self.line_color = config.line_color
        self._box_color = config.box_color
        self._box_border = config.box_border
        self._features_per_lane = 3
        self.features_per_lane = config.features_per_lane

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
