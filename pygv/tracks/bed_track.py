import os
import re
from collections import namedtuple
from typing import Any
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from .track import AnnotationTrack, Track


class _LaneRegistry(object):
    def __init__(self, offset=0, min_coord=None, max_coord=None, features=None):
        self.features = None
        if features is None:
            features = []
        self.offset = offset
        self.min_coord = min_coord
        self.max_coord = max_coord
        self.features = features


class BedTrack(AnnotationTrack):
    """
    If you're looking to visualize genomic features like genes and regulatory elements,
    you can utilize the :class:`~pygv.tracks.bed_track.BedTrack` in PyGV.
    If you desire greater control over the specific features to be plotted, such as filtering features by names,
    gene IDs, or transcript IDs, then the :class:`~pygv.tracks.gtf_track.GtfTrack` may be the preferred choice.

    Parameters
    ----------
    track : str
        Path to the input bed file. Index from tabix is optional, but
        when index presents, drawing will be much faster and less memory intensive.
        PyGV's BedTrack implementation supports the visualization of genomic features encoded in four BED variations:

        * Bed3: includes three columns: chromosome, start, and end.
        * Bed4: includes three columns: chromosome, start, end, and name.
        * Bed8: extends Bed3 with additional columns for name, score, strand, thickStart, and thickEnd.
        * Bed12: the most comprehensive, includes columns for itemRgb, blockCount, blockSizes, and blockStarts in addition to those found in Bed8.

    kwargs : Any
        height : float
            Height of each feature lane. If you have four feature lanes, and the height value
            is set as 0.25, then the final track will have identical overall height as other
            tracks. If you have multiple :class:`~pygv.tracks.bed_track.BedTrack` or other kind of :class:`~pygv.tracks.track.AnnotationTrack`, you can set
            the same height values to all these tracks, which makes sure all tracks have consistent
            appearances (which is also the default behavior).
        allowed_feature_lanes : Optional[int]
            See :attr:`~pygv.tracks.track.AnnotationTrack.allowed_feature_lanes`
        plot_thickness : bool
            See :attr:`~plot_thickness`
        padding_left : int or float
            See :attr:`~pygv.tracks.track.AnnotationTrack.padding_left`
        name : str
            Track name. :attr:`~pygv.tracks.track.Track.name`
        show_name : bool
            Show feature names. See :attr:`~pygv.tracks.track.AnnotationTrack.show_name`
        More kwargs can be seen here: :class:`~pygv.tracks.track.AnnotationTrack`

    Examples
    --------

    .. plot:: ../examples/plot_bed.py
    """

    def _get(self, chromosome, start, end):
        pass

    def _pysam_parser(self, chromosome, start, end):
        import pysam

        try:
            for row in self.bed_obj.fetch(
                chromosome, start, end, parser=pysam.asTuple()
            ):
                yield self._BedRecord._make(row)
        except ValueError:
            # in case no feature is available in that window
            return

    def _pd_parser(self, chromosome, start, end):
        for row in self.bed_obj.loc[
            np.logical_and(
                self.bed_obj.contig == chromosome,
                np.logical_and(self.bed_obj.start <= end, self.bed_obj.end >= start),
            ),
            :,
        ].values:
            yield self._BedRecord._make(row)

    def __init__(self, track, **kwargs: Any):
        if "height" not in kwargs:
            kwargs["height"] = 0.8
        super(BedTrack, self).__init__(track, **kwargs)
        is_bb = kwargs.pop("_is_bb", False)
        self._show_mode = ""
        self.show_mode = kwargs.pop("show_mode", "expanded")

        self._plot_thickness = 0
        self._plot_block = 0

        if not is_bb:
            if not os.path.exists(track):
                raise ValueError

            # parse bed file
            self.bed_file = track
            self.fields = (
                "contig",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts",
            )
            use_pysam = 1
            try:
                import pysam
            except ImportError:
                use_pysam = 0

            if use_pysam and track.endswith(".bed.gz"):
                if os.path.exists(track + ".tbi"):
                    use_pysam = 1
                else:
                    try:
                        pysam.tabix_index(track)
                        use_pysam = 1
                    except Exception as e:
                        use_pysam = 0
                        print(e)
            else:
                use_pysam = 0
            if use_pysam:
                self.bed_obj = pysam.TabixFile(track)
                tmp = pd.read_csv(track, sep="\t", header=None, comment="#", nrows=1)
                n_fields = tmp.shape[1]
                self._get = self._pysam_parser
            else:
                self.bed_obj = pd.read_csv(track, sep="\t", header=None, comment="#")
                n_fields = self.bed_obj.shape[1]
                self.bed_obj.columns = self.fields[:n_fields]
                self._get = self._pd_parser

            self._BedRecord = namedtuple("BedRecord", self.fields[:n_fields])

            if n_fields >= 8:
                self.plot_thickness = 1
                if "plot_thickness" in kwargs:
                    self.plot_thickness = kwargs.pop("plot_thickness")
            if n_fields == 12:
                self._plot_block = 1
        self._rgb_check = re.compile("(\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})")
        self.small_relative = 0
        self._block_line_width = 0
        self.block_line_height = kwargs.pop("block_line_height", 1)

        # override defaults
        if self.color is None:
            self.color = "#A1A1A1"
        if self.edge_color is None:
            self.edge_color = "#6E6E6E"

    @property
    def plot_thickness(self):
        """
        When valid starting (thickStart) and ending (thickEnd) positions are provided for a feature,
        as observed, for instance, in the region between the start and stop codon within gene displays,
        this specific region will be illustrated with a thicker line. In situations where there is no
        valid `thick` section, whether due to the absence of `thickStart` and `thickEnd` positions
        (as in the case of providing a bed3 file) or when the length of the feature is less than 1,
        the :code:`plot_thickness` property will be automatically set to :code:`False`.
        """
        return self._plot_thickness

    @plot_thickness.setter
    def plot_thickness(self, value):
        try:
            self._plot_thickness = bool(value)
        except Exception as e:
            print(e)

    @property
    def show_mode(self):
        """
        Collapse overlapping features (collapsed) or keep them separately (expanded) for plotting.
        """
        return self._show_mode

    @show_mode.setter
    def show_mode(self, value):
        if value not in ["expanded", "collapsed"]:
            raise ValueError("show_mode must be either expanded or collapsed")
        else:
            self._show_mode = value

    @property
    def block_line_width(self):
        """
        Line/edge width for blocks.
        """
        return self._block_line_width

    @block_line_width.setter
    def block_line_width(self, value):
        try:
            self._block_line_width = float(value)
        except Exception as e:
            raise e

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
        region_len = end - start
        text_padding = (
            region_len * self.padding_left
            if 0 < self.padding_left < 1
            else self.padding_left
        )
        self._lane_registries = []
        added = set()
        for interval in self._get(chromosome=chromosome, start=start, end=end):
            active_lane = None
            start_loc = int(interval.start)
            end_loc = int(interval.end)
            visible_start = max(start_loc, start)
            visible_end = min(end_loc, end)

            if self._hide_visual_dup:
                k = (visible_start, visible_end, interval.strand)
                if k in added:
                    continue
                added.add(k)

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
                    if (
                        lr.max_coord < start_loc - text_padding
                        or self.show_mode == "collapsed"
                    ):
                        active_lane = lr.offset
                        lr.min_coord = min(lr.min_coord, start_loc)
                        lr.max_coord = max(lr.max_coord, end_loc)
                        lr.features.append(interval)
                        break
                    else:
                        active_lane = None

            if (
                type(self.allowed_feature_lanes) is int
                and len(self._lane_registries) >= self.allowed_feature_lanes
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
        if len(self._lane_registries) == 1:
            self._lane_registries.append(
                _LaneRegistry(
                    offset=len(self._lane_registries),
                    min_coord=start,
                    max_coord=end,
                    features=[],
                )
            )

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw track

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the region of interest/window, 0-based
        end : int
            end of the region of interest/window, 0-based
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(BedTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        import matplotlib.pyplot as plt

        fig = plt.gcf()
        self._ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        self._ax.set_xlim((start, end))

        self._small_relative = 0.004 * (end - start)
        for lane in self._lane_registries:
            empty_lane = True
            active_lane = lane.offset
            real_active_line = (self._patch_height + self._lane_space) * active_lane
            for interval in lane.features:
                color = self.color
                empty_lane = False
                start_loc = int(interval.start)
                end_loc = int(interval.end)
                visible_start = max(start_loc, start)
                visible_end = min(end_loc, end)
                plot_thickness = self.plot_thickness

                if "itemRgb" in dir(interval):
                    try:
                        m = self._rgb_check.match(interval.itemRgb)
                    except:
                        m = False
                    if m:
                        r, g, b = map(int, m.groups())
                        color = colors.to_hex((r / 255, g / 255, b / 255))

                if not self._plot_block:
                    # first, draw an invisible line for the determination of correct ylim
                    self._ax.plot(
                        (
                            start_loc if start_loc >= start else visible_start,
                            end_loc - 1 if end_loc <= end else visible_end,
                        ),
                        (-1 * real_active_line, -1 * real_active_line),
                        color=self._line_color,
                        linewidth=self.line_width,
                        alpha=0,
                        clip_on=True,
                        zorder=-1,
                    )
                    # then draw the visible box
                    rec = Rectangle(
                        xy=(
                            start_loc,
                            -1 * real_active_line - (self._patch_height / 2),
                        ),
                        width=end_loc - start_loc,
                        height=self._patch_height,
                        edgecolor=self.edge_color,
                        facecolor=color,
                        linewidth=self.block_line_height,
                    )
                    self._ax.add_patch(rec)
                else:
                    # init
                    exon_starts = list(
                        map(int, filter(None, interval.blockStarts.split(",")))
                    )
                    exon_sizes = list(
                        map(int, filter(None, interval.blockSizes.split(",")))
                    )
                    # if the record has matched number of exon starts and sizes, and this number is
                    # consistent with the blockCount, then we try to plot each block/exon
                    if len(exon_starts) == len(exon_sizes) == int(interval.blockCount):
                        self._ax.plot(
                            (
                                start_loc if start_loc >= start else visible_start,
                                end_loc - 1 if end_loc <= end else visible_end,
                            ),
                            (-1 * real_active_line, -1 * real_active_line),
                            color=self._line_color,
                            linewidth=self.line_width,
                            clip_on=True,
                        )

                        # plot small arrows over the backbone
                        if end_loc - start_loc > self._small_relative:
                            pos = np.arange(
                                visible_start + self._small_relative,
                                visible_end + self._small_relative,
                                int(self._arrow_interval * self._small_relative),
                            )
                            for xpos in pos:
                                self._plot_gene_direction(
                                    ax, xpos, -1 * real_active_line, interval.strand
                                )

                        patches = []
                        try:
                            thick_start = int(interval.thickStart)
                            thick_end = int(interval.thickEnd)
                            if thick_end - thick_start == 0:
                                plot_thickness = 0
                        except:
                            thick_start = None
                            thick_end = None
                            plot_thickness = 0
                        for i, (s, size) in enumerate(zip(exon_starts, exon_sizes)):
                            x = s + start_loc
                            if x < end or (x + size) > start:
                                # in case of overflow
                                adjusted_size = size
                                plot_end = x + adjusted_size
                                if plot_end > visible_end:
                                    adjusted_size -= plot_end - visible_end
                                if plot_thickness:
                                    if x < thick_start:
                                        begin_of_thickness = min(
                                            adjusted_size, max(thick_start - x, 0)
                                        )
                                        # thinner part
                                        p = Rectangle(
                                            xy=(
                                                x,
                                                -1 * real_active_line
                                                - (self._patch_height / 4),
                                            ),
                                            width=begin_of_thickness,
                                            clip_on=True,
                                            height=self._patch_height / 2,
                                        )
                                        patches.append(p)
                                        # the thick part is larger than the current block
                                        if x + size < thick_end:
                                            # thicker part
                                            remaining_length = abs(
                                                adjusted_size - begin_of_thickness
                                            )
                                            if remaining_length > 0:
                                                p = Rectangle(
                                                    xy=(
                                                        x + begin_of_thickness,
                                                        -1 * real_active_line
                                                        - (self._patch_height / 2),
                                                    ),
                                                    width=remaining_length,
                                                    clip_on=True,
                                                    height=self._patch_height,
                                                )
                                                patches.append(p)
                                        else:
                                            p = Rectangle(
                                                xy=(
                                                    x + begin_of_thickness,
                                                    -1 * real_active_line
                                                    - (self._patch_height / 2),
                                                ),
                                                width=max(
                                                    thick_end - x - begin_of_thickness,
                                                    0,
                                                ),
                                                clip_on=True,
                                                height=self._patch_height,
                                            )
                                            patches.append(p)
                                            # remaining part
                                            p = Rectangle(
                                                xy=(
                                                    thick_end,
                                                    -1 * real_active_line
                                                    - (self._patch_height / 4),
                                                ),
                                                width=max(x + size - thick_end, 0),
                                                clip_on=True,
                                                height=self._patch_height / 2,
                                            )
                                            patches.append(p)
                                    elif x < thick_end:
                                        if x + size < thick_end:
                                            p = Rectangle(
                                                xy=(
                                                    x,
                                                    -1 * real_active_line
                                                    - (self._patch_height / 2),
                                                ),
                                                width=adjusted_size,
                                                clip_on=True,
                                                height=self._patch_height,
                                            )
                                            patches.append(p)
                                        else:
                                            # thinner part
                                            end_of_thickness = thick_end
                                            if end_of_thickness < visible_end:
                                                p = Rectangle(
                                                    xy=(
                                                        end_of_thickness,
                                                        -1 * real_active_line
                                                        - (self._patch_height / 4),
                                                    ),
                                                    width=max(
                                                        size - end_of_thickness + x, 0
                                                    ),
                                                    clip_on=True,
                                                    height=self._patch_height / 2,
                                                )
                                                # thicker part
                                                patches.append(p)
                                            remaining_length = abs(end_of_thickness - x)
                                            if remaining_length > 0:
                                                p = Rectangle(
                                                    xy=(
                                                        x,
                                                        -1 * real_active_line
                                                        - (self._patch_height / 2),
                                                    ),
                                                    width=remaining_length,
                                                    clip_on=True,
                                                    height=self._patch_height,
                                                )
                                                patches.append(p)
                                    else:
                                        p = Rectangle(
                                            xy=(
                                                x,
                                                -1 * real_active_line
                                                - (self._patch_height / 4),
                                            ),
                                            width=adjusted_size,
                                            clip_on=True,
                                            height=self._patch_height / 2,
                                        )
                                        patches.append(p)
                                else:
                                    p = Rectangle(
                                        xy=(
                                            x,
                                            -1 * real_active_line
                                            - (self._patch_height / 2),
                                        ),
                                        width=adjusted_size,
                                        clip_on=True,
                                        height=self._patch_height,
                                    )
                                    patches.append(p)

                        self._ax.add_collection(
                            PatchCollection(
                                patches,
                                edgecolors=self.edge_color,
                                facecolors=color,
                                linewidths=self.block_line_height,
                                zorder=100,
                                clip_on=True,
                            )
                        )
                    else:  # otherwise, we plot a single bar
                        rec = Rectangle(
                            xy=(
                                start_loc,
                                -1 * real_active_line - self._patch_height / 2,
                            ),
                            width=end_loc - start_loc,
                            height=self._patch_height,
                            edgecolor=self.edge_color,
                            facecolor=color,
                            linewidth=self.block_line_height,
                            **kwargs,
                        )
                        self._ax.add_patch(rec)

                if "name" in dir(interval) and self.show_name:
                    if (
                        start_loc > start
                        and "strand" in dir(interval)
                        and interval.strand == "+"
                    ):
                        self._ax.text(
                            x=start_loc - self._small_relative,
                            y=-1 * real_active_line,
                            color=self.font_color,
                            size=self.font_size,
                            s=interval.name,
                            ha="right",
                            va="center",
                            clip_on=True,
                            zorder=101,
                        )
                    elif (
                        end_loc < end
                        and "strand" in dir(interval)
                        and interval.strand == "-"
                    ):
                        self._ax.text(
                            x=end_loc + self._small_relative,
                            y=-1 * real_active_line,
                            color=self.font_color,
                            size=self.font_size,
                            s=interval.name,
                            ha="left",
                            va="center",
                            clip_on=True,
                            zorder=101,
                        )
                    else:
                        self._ax.text(
                            x=(visible_end + visible_start) / 2,
                            y=-1 * real_active_line,
                            color=self.font_color,
                            size=self.font_size,
                            s=interval.name,
                            ha="center",
                            va="center",
                            clip_on=True,
                            bbox=dict(
                                boxstyle="round",
                                fc="w",
                                alpha=self.font_box_alpha,
                                lw=0.1,
                            ),
                            zorder=101,
                        )

            if empty_lane:
                # first, draw an invisible line for the determination of correct ylim
                self._ax.plot(
                    (start, start + 1),
                    (-1 * real_active_line, -1 * real_active_line),
                    color=self._line_color,
                    linewidth=self.line_width,
                    alpha=0,
                    clip_on=True,
                    zorder=-1,
                )
        self._ax.set_yticks([])
        # remove minor ticks
        self._ax.set_yticks([], minor=True)

        if index != 0:
            # remove major ticks
            self._ax.set_xticks([])
            # remove minor ticks
            self._ax.set_xticks([], minor=True)
            # self.ax.margins(0)


class BedPETrack(AnnotationTrack):
    """
    If you're looking to visualize genomic interactions like proximity captured by Hi-C,
    you can utilize the :class:`~pygv.tracks.bed_track.BedPETrack` in PyGV.

    Parameters
    ----------
    track : str
        Path to the input BEDPE file. Index from tabix is optional, but
        when index presents, drawing will be much faster and less memory intensive.
        Only the first six columns in the file will be used
        (chromosome, start, and end for the two anchors)
    kwargs : dict
        name : str
            Track name. :attr:`~pygv.tracks.track.Track.name`
        show_name : bool
            Show feature names. See :attr:`~pygv.tracks.track.AnnotationTrack.show_name`
        flip : bool
            Flip the arcs vertically
        More kwargs can be seen here: :class:`~pygv.tracks.track.AnnotationTrack`

    Examples
    --------

    .. plot:: ../examples/plot_bedpe.py
    """

    def _get(self, chromosome, start, end):
        pass

    def _pysam_parser(self, chromosome, start, end):
        import pysam

        try:
            for row in self.bed_obj.fetch(
                chromosome, start, end, parser=pysam.asTuple()
            ):
                yield self._BedPERecord._make(row[: len(self._BedPERecord._fields)])
        except ValueError:
            # in case no feature is available in that window
            return

    def _pd_parser(self, chromosome, start, end):
        for row in self.bed_obj.loc[
            np.logical_and(
                self.bed_obj.chrom1 == chromosome,
                np.logical_and(self.bed_obj.start1 <= end, self.bed_obj.end1 >= start),
            ),
            :,
        ].values:
            yield self._BedPERecord._make(row[: len(self._BedPERecord._fields)])

    def __init__(self, track, **kwargs: Any):
        super(BedPETrack, self).__init__(track, **kwargs)
        if not os.path.exists(track):
            raise ValueError

        # parse bed file
        self.bed_file = track
        self.fields = (
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "name",
            "score",
            "strand1",
            "strand2",
            "others",
        )
        use_pysam = 1
        try:
            import pysam
        except ImportError:
            use_pysam = 0

        if use_pysam and track.endswith(".bedpe.gz"):
            if os.path.exists(track + ".tbi"):
                use_pysam = 1
            else:
                try:
                    pysam.tabix_index(track, preset="bed")
                    use_pysam = 1
                except Exception as e:
                    use_pysam = 0
                    print(e)
        if use_pysam:
            self.bed_obj = pysam.TabixFile(track)
            tmp = pd.read_csv(track, sep="\t", header=None, comment="#", nrows=1)
            n_fields = tmp.shape[1]
            self._get = self._pysam_parser
        else:
            self.bed_obj = pd.read_csv(track, sep="\t", header=None, comment="#")
            n_fields = min(self.bed_obj.shape[1], len(self.fields))
            self.bed_obj.columns = self.fields[:n_fields]
            self._get = self._pd_parser

        self._BedPERecord = namedtuple("BedPERecord", self.fields[:n_fields])
        self.small_relative = 0

        # override defaults
        if self.color is None:
            self.color = "#A1A1A1"
        if self.edge_color is None:
            self.edge_color = "#6E6E6E"

        self.flip_arc = kwargs.get("flip", False)

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
        self._lane_registries = [[]]
        added = set()
        for interval in self._get(chromosome=chromosome, start=start, end=end):
            start_loc = min(int(interval.start1), int(interval.start2))
            end_loc = max(int(interval.end1), int(interval.end2))
            visible_start = max(start_loc, start)
            visible_end = min(end_loc, end)

            if self._hide_visual_dup:
                k = (visible_start, visible_end, interval.strand)
                if k in added:
                    continue
                added.add(k)
            self._lane_registries[0].append(
                (
                    int(interval.start1),
                    int(interval.end1),
                    int(interval.start2),
                    int(interval.end2),
                )
            )

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw track

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the region of interest/window, 0-based
        end : int
            end of the region of interest/window, 0-based
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), :code:`index==0`, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(BedPETrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        from matplotlib.patches import Arc

        self._ax.set_xlim((start, end))

        self._small_relative = 0.004 * (end - start)
        for pair in self._lane_registries[0]:
            a1_start_loc, a1_end_loc, a2_start_loc, a2_end_loc = pair
            a2_mid = (a2_end_loc + a2_start_loc) / 2
            a1_mid = (a1_end_loc + a1_start_loc) / 2
            self._ax.plot(
                (a1_start_loc, a1_end_loc),
                (0, 0),
                color=self.color,
                lw=self.line_width,
                clip_on=False,
            )
            self._ax.plot(
                (a2_start_loc, a2_end_loc), (0, 0), color=self.color, lw=self.line_width
            )
            if a1_start_loc < a2_start_loc:
                arc_length = a2_mid - a1_mid
                x = a1_mid + arc_length / 2
            else:
                arc_length = a1_mid - a2_mid
                x = a2_mid + arc_length / 2
            arc = Arc(
                (x, 0),
                arc_length,
                arc_length / 2,
                theta2=180,
                color=self.color,
                alpha=self.alpha,
            )
            ax.add_patch(arc)

        self._ax.set_yticks([])
        # remove minor ticks
        self._ax.set_yticks([], minor=True)
        self._ax.set_ylim((0, None))
        self._ax.autoscale_view()
        if index != 0:
            # remove major ticks
            self._ax.set_xticks([])
            # remove minor ticks
            self._ax.set_xticks([], minor=True)

        if self.flip_arc:
            self._ax.invert_yaxis()


class ConnectionArcTrack(BedTrack):
    """
    Similar to :class:`~pygv.tracks.bed_track.BedPETrack`, :class:`~pygv.tracks.bed_track.ConnectionArcTrack`
    can be used to visualize genomic interactions like enhancer-promoter interaction.
    The main difference between these two tracks is that :class:`~pygv.tracks.bed_track.ConnectionArcTrack`
    draws a directed arrow from the source to the target.

    Parameters
    ----------
    track : str
        Path to the input BEDPE file. Index from tabix is optional, but
        when index presents, drawing will be much faster and less memory intensive.
        Only the first six columns in the file will be used
        (chromosome, start, and end for the two anchors)
    kwargs : dict
        name : str
            Track name. :attr:`~pygv.tracks.track.Track.name`
        show_name : bool
            Show feature names. See :attr:`~pygv.tracks.track.AnnotationTrack.show_name`
        flip : bool
            Flip the arcs vertically
        More kwargs can be seen here: :class:`~pygv.tracks.track.AnnotationTrack`

    Examples
    --------

    .. plot:: ../examples/plot_bedpe.py
    """

    def __init__(self, track, **kwargs):
        super(ConnectionArcTrack, self).__init__(track=track, **kwargs)
        if self.color is None:
            self.color = "#FFD900"
        if self.edge_color is None:
            self.edge_color = "#FFBD00"

        if not os.path.exists(track):
            raise IOError("File not found: {}".format(track))

        self.flip_arc = kwargs.get("flip", False)

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
        pass

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """

        Parameters
        ----------
        chromosome :
        start :
        end :
        ax :
        index :
        kwargs :

        Returns
        -------

        """
        super(ConnectionArcTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        arrow_style = kwargs.pop("arrow_style", "->")
        connection_style = kwargs.pop("connection_style", None)
        rad = kwargs.pop("rad", 0.2)
        starts = []
        ends = []
        for interval in self._get(chromosome=chromosome, start=start, end=end):
            anchor_start = int(interval.start)
            anchor_end = int(interval.end)
            anchor_mid = (anchor_start + anchor_end) / 2
            # ugly naming...
            target_start = int(interval.score)
            target_end = int(interval.strand)

            # check whether the target is at least partially in the visible window
            if start <= target_end <= end or start <= target_start <= end:
                if anchor_mid <= target_start:
                    if target_start >= start:
                        if connection_style is None:
                            connection_style = "arc3,rad=-{rad}".format(rad=rad)
                        ax.annotate(
                            "",
                            xy=(target_start, 0),
                            xycoords="data",
                            xytext=(anchor_mid, self._height),
                            textcoords="data",
                            arrowprops=dict(
                                arrowstyle=arrow_style,
                                connectionstyle=connection_style,
                                linewidth=self.line_width,
                            ),
                        )
                    else:
                        continue
                else:
                    if target_end <= end:
                        if connection_style is None:
                            connection_style = "arc3,rad={rad}".format(rad=rad)
                        ax.annotate(
                            "",
                            xy=(target_end, 0),
                            xycoords="data",
                            xytext=(anchor_mid, self._height),
                            textcoords="data",
                            arrowprops=dict(
                                arrowstyle=arrow_style,
                                connectionstyle=connection_style,
                                linewidth=self.line_width,
                            ),
                        )
                    else:
                        continue
                p = Rectangle(
                    xy=(anchor_start, 0),
                    width=anchor_end - anchor_start,
                    edgecolor=self.edge_color,
                    clip_on=False,
                    height=self._height,
                    facecolor=self.color,
                    alpha=self._alpha,
                    linewidth=self.line_width,
                    zorder=100,
                )
                self._ax.add_patch(p)
                starts.append(min(anchor_mid, target_start, target_end))
                ends.append(max(anchor_mid, target_start, target_end))

        self._ax.yaxis.set_ticks([])
        self._ax.set_xlim((start, end))
        try:
            self._ax.ticklabel_format(style="plain", useOffset=False)
        except:
            pass
        # self.ax.autoscale_view()
        self._ax.set_ylim((-0.05, 3))

        if self.flip_arc:
            self._ax.invert_yaxis()
