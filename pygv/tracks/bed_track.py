import os
import re
from collections import namedtuple
from typing import Any, Optional

import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from pydantic import AliasChoices, Field, PrivateAttr

from .track import AnnotationTrack
from .types import Color, ShowMode


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
    """Visualize genomic features from BED files (BED3/4/8/12)."""

    height: float = Field(
        default=0.8,
        gt=0,
        kw_only=True,
        description=(
            "Height of each feature lane. If you have four feature lanes and height "
            "is 0.25, the final track has the same overall height as a unit-height track."
        ),
    )
    show_mode: ShowMode = Field(
        default="expanded",
        kw_only=True,
        description=(
            "Collapse overlapping features (`collapsed`) or keep them separately "
            "(`expanded`) for plotting."
        ),
    )
    plot_thickness: Optional[bool] = Field(
        default=None,
        kw_only=True,
        description=(
            "When thickStart/thickEnd are present (e.g. CDS), draw that region with "
            "a thicker box. Automatically enabled for BED8+ unless set explicitly."
        ),
    )
    block_line_height: float = Field(
        default=1,
        kw_only=True,
        description="Line/edge width for blocks",
    )

    _bed_file: Optional[str] = PrivateAttr(default=None)
    _fields: tuple = PrivateAttr(default=())
    _bed_obj: Any = PrivateAttr(default=None)
    _BedRecord: Any = PrivateAttr(default=None)
    _parser: Any = PrivateAttr(default=None)
    _rgb_check: Any = PrivateAttr(default=None)

    def _get(self, chromosome, start, end):
        if self._parser is None:
            return
        yield from self._parser(chromosome, start, end)

    def _pysam_parser(self, chromosome, start, end):
        import pysam

        try:
            for row in self._bed_obj.fetch(
                chromosome, start, end, parser=pysam.asTuple()
            ):
                yield self._BedRecord._make(row)
        except ValueError:
            # in case no feature is available in that window
            return

    def _pd_parser(self, chromosome, start, end):
        for row in self._bed_obj.loc[
            np.logical_and(
                self._bed_obj.contig == chromosome,
                np.logical_and(self._bed_obj.start <= end, self._bed_obj.end >= start),
            ),
            :,
        ].values:
            yield self._BedRecord._make(row)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        self._open_source()

    def _open_source(self) -> None:
        if not os.path.exists(self.track):
            raise ValueError

        self._bed_file = self.track
        self._fields = (
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

        if use_pysam and self.track.endswith(".bed.gz"):
            if os.path.exists(self.track + ".tbi"):
                use_pysam = 1
            else:
                try:
                    pysam.tabix_index(self.track)
                    use_pysam = 1
                except Exception as e:
                    use_pysam = 0
                    print(e)
        else:
            use_pysam = 0
        if use_pysam:
            self._bed_obj = pysam.TabixFile(self.track)
            tmp = pd.read_csv(self.track, sep="\t", header=None, comment="#", nrows=1)
            n_fields = tmp.shape[1]
            self._parser = self._pysam_parser
        else:
            self._bed_obj = pd.read_csv(
                self.track, sep="\t", header=None, comment="#"
            )
            n_fields = self._bed_obj.shape[1]
            self._bed_obj.columns = self._fields[:n_fields]
            self._parser = self._pd_parser

        self._BedRecord = namedtuple("BedRecord", self._fields[:n_fields])

        if n_fields >= 8:
            if self.plot_thickness is None:
                self.plot_thickness = True
        elif self.plot_thickness is None:
            self.plot_thickness = False
        if n_fields == 12:
            self._plot_block = 1
        self._rgb_check = re.compile(r"(\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})")
        self._small_relative = 0

        if self.color is None:
            self.color = "#A1A1A1"
        if self.edge_color is None:
            self.edge_color = "#6E6E6E"

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

            if self.hide_visual_dup:
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
            real_active_line = (self.patch_height + self.lane_space) * active_lane
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
                        color=self.line_color,
                        linewidth=self.line_width,
                        alpha=0,
                        clip_on=True,
                        zorder=-1,
                    )
                    # then draw the visible box
                    rec = Rectangle(
                        xy=(
                            start_loc,
                            -1 * real_active_line - (self.patch_height / 2),
                        ),
                        width=end_loc - start_loc,
                        height=self.patch_height,
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
                            color=self.line_color,
                            linewidth=self.line_width,
                            clip_on=True,
                        )

                        # plot small arrows over the backbone
                        if end_loc - start_loc > self._small_relative:
                            pos = np.arange(
                                visible_start + self._small_relative,
                                visible_end + self._small_relative,
                                int(self.arrow_interval * self._small_relative),
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
                                                - (self.patch_height / 4),
                                            ),
                                            width=begin_of_thickness,
                                            clip_on=True,
                                            height=self.patch_height / 2,
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
                                                        - (self.patch_height / 2),
                                                    ),
                                                    width=remaining_length,
                                                    clip_on=True,
                                                    height=self.patch_height,
                                                )
                                                patches.append(p)
                                        else:
                                            p = Rectangle(
                                                xy=(
                                                    x + begin_of_thickness,
                                                    -1 * real_active_line
                                                    - (self.patch_height / 2),
                                                ),
                                                width=max(
                                                    thick_end - x - begin_of_thickness,
                                                    0,
                                                ),
                                                clip_on=True,
                                                height=self.patch_height,
                                            )
                                            patches.append(p)
                                            # remaining part
                                            p = Rectangle(
                                                xy=(
                                                    thick_end,
                                                    -1 * real_active_line
                                                    - (self.patch_height / 4),
                                                ),
                                                width=max(x + size - thick_end, 0),
                                                clip_on=True,
                                                height=self.patch_height / 2,
                                            )
                                            patches.append(p)
                                    elif x < thick_end:
                                        if x + size < thick_end:
                                            p = Rectangle(
                                                xy=(
                                                    x,
                                                    -1 * real_active_line
                                                    - (self.patch_height / 2),
                                                ),
                                                width=adjusted_size,
                                                clip_on=True,
                                                height=self.patch_height,
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
                                                        - (self.patch_height / 4),
                                                    ),
                                                    width=max(
                                                        size - end_of_thickness + x, 0
                                                    ),
                                                    clip_on=True,
                                                    height=self.patch_height / 2,
                                                )
                                                # thicker part
                                                patches.append(p)
                                            remaining_length = abs(end_of_thickness - x)
                                            if remaining_length > 0:
                                                p = Rectangle(
                                                    xy=(
                                                        x,
                                                        -1 * real_active_line
                                                        - (self.patch_height / 2),
                                                    ),
                                                    width=remaining_length,
                                                    clip_on=True,
                                                    height=self.patch_height,
                                                )
                                                patches.append(p)
                                    else:
                                        p = Rectangle(
                                            xy=(
                                                x,
                                                -1 * real_active_line
                                                - (self.patch_height / 4),
                                            ),
                                            width=adjusted_size,
                                            clip_on=True,
                                            height=self.patch_height / 2,
                                        )
                                        patches.append(p)
                                else:
                                    p = Rectangle(
                                        xy=(
                                            x,
                                            -1 * real_active_line
                                            - (self.patch_height / 2),
                                        ),
                                        width=adjusted_size,
                                        clip_on=True,
                                        height=self.patch_height,
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
                                -1 * real_active_line - self.patch_height / 2,
                            ),
                            width=end_loc - start_loc,
                            height=self.patch_height,
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
                    color=self.line_color,
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
    """Visualize genomic interactions from BEDPE files."""

    flip_arc: bool = Field(
        default=False,
        kw_only=True,
        validation_alias=AliasChoices("flip_arc", "flip"),
        description="Flip the arcs vertically",
    )
    highlight_link_color: Color = Field(
        default="#D62728",
        kw_only=True,
        description="Color used for highlighted links",
    )
    highlight_link_alpha: float = Field(
        default=1.0,
        ge=0,
        le=1,
        kw_only=True,
        description="Alpha used for highlighted links",
    )
    highlight_link_line_width: Optional[float] = Field(
        default=None,
        kw_only=True,
        description="Line width used for highlighted links",
    )
    highlight_links: Any = Field(
        default=(),
        kw_only=True,
        description=(
            "Links to highlight. Each item can be a BEDPE name string or a "
            "coordinate tuple (chromosome, start1, end1, start2, end2)."
        ),
    )

    _bed_file: Optional[str] = PrivateAttr(default=None)
    _fields: tuple = PrivateAttr(default=())
    _bed_obj: Any = PrivateAttr(default=None)
    _BedPERecord: Any = PrivateAttr(default=None)
    _parser: Any = PrivateAttr(default=None)
    _highlight_links: set = PrivateAttr(default_factory=set)

    def _get(self, chromosome, start, end):
        if self._parser is None:
            return
        yield from self._parser(chromosome, start, end)

    def _pysam_parser(self, chromosome, start, end):
        import pysam

        try:
            for row in self._bed_obj.fetch(
                chromosome, start, end, parser=pysam.asTuple()
            ):
                yield self._BedPERecord._make(row[: len(self._BedPERecord._fields)])
        except ValueError:
            # in case no feature is available in that window
            return

    def _pd_parser(self, chromosome, start, end):
        for row in self._bed_obj.loc[
            np.logical_and(
                self._bed_obj.chrom1 == chromosome,
                np.logical_and(self._bed_obj.start1 <= end, self._bed_obj.end1 >= start),
            ),
            :,
        ].values:
            yield self._BedPERecord._make(row[: len(self._BedPERecord._fields)])

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        if not os.path.exists(self.track):
            raise ValueError

        self._bed_file = self.track
        self._fields = (
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

        if use_pysam and self.track.endswith(".bedpe.gz"):
            if os.path.exists(self.track + ".tbi"):
                use_pysam = 1
            else:
                try:
                    pysam.tabix_index(self.track, preset="bed")
                    use_pysam = 1
                except Exception as e:
                    use_pysam = 0
                    print(e)
        else:
            use_pysam = 0
        if use_pysam:
            self._bed_obj = pysam.TabixFile(self.track)
            tmp = pd.read_csv(self.track, sep="\t", header=None, comment="#", nrows=1)
            n_fields = tmp.shape[1]
            self._parser = self._pysam_parser
        else:
            self._bed_obj = pd.read_csv(
                self.track, sep="\t", header=None, comment="#"
            )
            n_fields = min(self._bed_obj.shape[1], len(self._fields))
            self._bed_obj.columns = self._fields[:n_fields]
            self._parser = self._pd_parser

        self._BedPERecord = namedtuple("BedPERecord", self._fields[:n_fields])
        self._small_relative = 0

        if self.color is None:
            self.color = "#A1A1A1"
        if self.edge_color is None:
            self.edge_color = "#6E6E6E"
        if self.highlight_link_line_width is None:
            self.highlight_link_line_width = self.line_width * 1.5
        self.set_highlight_links(self.highlight_links)

    @staticmethod
    def _normalize_link_coords(start1, end1, start2, end2, chromosome=None):
        first_anchor = (int(start1), int(end1))
        second_anchor = (int(start2), int(end2))
        left_anchor, right_anchor = sorted((first_anchor, second_anchor))
        normalized_chromosome = None if chromosome is None else str(chromosome)
        return (
            "coords",
            normalized_chromosome,
            left_anchor[0],
            left_anchor[1],
            right_anchor[0],
            right_anchor[1],
        )

    def _normalize_highlight_link(self, link):
        if isinstance(link, str):
            return ("name", str(link))
        if isinstance(link, (tuple, list)) and len(link) == 5:
            chromosome, start1, end1, start2, end2 = link
            try:
                return self._normalize_link_coords(
                    start1, end1, start2, end2, chromosome=chromosome
                )
            except (TypeError, ValueError):
                return None
        if isinstance(link, (tuple, list)) and len(link) == 4:
            try:
                return self._normalize_link_coords(*link)
            except (TypeError, ValueError):
                return None
        return None

    def set_highlight_links(self, links):
        """
        Set highlighted links for this BEDPE track.

        Parameters
        ----------
        links : str or iterable
            A BEDPE `name`, or an iterable of names and/or coordinate tuples
            `(chromosome, start1, end1, start2, end2)`.
        """
        if links is None:
            self._highlight_links = set()
            return

        if isinstance(links, str):
            links = (links,)

        normalized_links = set()
        for link in links:
            normalized = self._normalize_highlight_link(link)
            if normalized is not None:
                normalized_links.add(normalized)
        self._highlight_links = normalized_links

    def add_highlight_link(self, link):
        """
        Add one highlighted link by BEDPE name or anchor coordinates.
        """
        normalized = self._normalize_highlight_link(link)
        if normalized is not None:
            self._highlight_links.add(normalized)

    def clear_highlight_links(self):
        """
        Remove all highlighted links.
        """
        self._highlight_links = set()

    def _is_highlight_link(self, chromosome, name, start1, end1, start2, end2):
        if len(self._highlight_links) == 0:
            return False
        if name is not None and ("name", str(name)) in self._highlight_links:
            return True
        return self._normalize_link_coords(
            start1, end1, start2, end2, chromosome=chromosome
        ) in self._highlight_links or self._normalize_link_coords(
            start1, end1, start2, end2
        ) in self._highlight_links

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

            if self.hide_visual_dup:
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
                    getattr(interval, "name", None),
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
            if len(pair) >= 5:
                a1_start_loc, a1_end_loc, a2_start_loc, a2_end_loc, pair_name = pair
            else:
                a1_start_loc, a1_end_loc, a2_start_loc, a2_end_loc = pair
                pair_name = None

            is_highlight = self._is_highlight_link(
                chromosome,
                pair_name,
                a1_start_loc,
                a1_end_loc,
                a2_start_loc,
                a2_end_loc,
            )
            plot_color = self.highlight_link_color if is_highlight else self.color
            plot_alpha = self.highlight_link_alpha if is_highlight else self.alpha
            plot_line_width = (
                self.highlight_link_line_width if is_highlight else self.line_width
            )
            a2_mid = (a2_end_loc + a2_start_loc) / 2
            a1_mid = (a1_end_loc + a1_start_loc) / 2
            self._ax.plot(
                (a1_start_loc, a1_end_loc),
                (0, 0),
                color=plot_color,
                lw=plot_line_width,
                alpha=plot_alpha,
                clip_on=False,
            )
            self._ax.plot(
                (a2_start_loc, a2_end_loc),
                (0, 0),
                color=plot_color,
                lw=plot_line_width,
                alpha=plot_alpha,
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
                color=plot_color,
                alpha=plot_alpha,
                linewidth=plot_line_width,
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
    """Directed arcs from source to target (BED-like input)."""

    color: Color = Field(default="#FFD900", kw_only=True)
    edge_color: Color = Field(default="#FFBD00", kw_only=True)
    flip_arc: bool = Field(
        default=False,
        kw_only=True,
        validation_alias=AliasChoices("flip_arc", "flip"),
        description="Flip the arcs vertically",
    )
    arrow_style: str = Field(
        default="->", kw_only=True, description="Matplotlib arrow style"
    )
    connection_style: Optional[str] = Field(
        default=None,
        kw_only=True,
        description="Matplotlib connection style; default is an arc using `rad`",
    )
    rad: float = Field(
        default=0.2, kw_only=True, description="Arc radius used when connection_style is None"
    )

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
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
        arrow_style = self.arrow_style
        connection_style = self.connection_style
        rad = self.rad
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
                            xytext=(anchor_mid, self.height),
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
                            xytext=(anchor_mid, self.height),
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
                    height=self.height,
                    facecolor=self.color,
                    alpha=self.alpha,
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
