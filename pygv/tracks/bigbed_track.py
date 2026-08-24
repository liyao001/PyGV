import os
import re
import warnings
from collections import namedtuple
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig
from matplotlib.lines import Line2D
from pydantic import Field, PrivateAttr

from .bed_track import BedTrack
from .track import NumericalTrack
from .types import Color


class UCSCMutationTrack(NumericalTrack):
    """Lollipop plot from UCSC-style mutational bigBed files."""

    track: str = Field(description="Path or URL to a bigBed file")
    line_color: Color = Field(default="red", kw_only=True, description="Stem color")
    apply_color_gradient: bool = Field(
        default=False,
        kw_only=True,
        description="Apply a color gradient to mutation markers",
    )
    color_map: Any = Field(
        default_factory=lambda: plt.cm.Reds,
        kw_only=True,
        description="Colormap for markers when apply_color_gradient is True",
    )
    normalizer: Any = Field(
        default_factory=lambda: matplotlib.colors.Normalize,
        kw_only=True,
        description="Matplotlib normalizer class for the color gradient",
    )

    _bb: Any = PrivateAttr(default=None)
    _filters: dict = PrivateAttr(default_factory=dict)
    _filter_supported_fields: set = PrivateAttr(default_factory=lambda: {"MAF", "ID"})
    _color_map: Any = PrivateAttr(default=None)

    def __init__(self, track: str, **data: Any) -> None:
        super().__init__(track=track, **data)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        if not os.path.exists(self.track) and not self.track.startswith("http"):
            raise ValueError
        try:
            self._color_map = matplotlib.colors.ListedColormap(
                self.color_map(np.linspace(0, 1, 20))[:-5, :-1]
            )
        except Exception:
            self._color_map = self.color_map
        self._bb = pyBigWig.open(self.track)
        if not self._bb.isBigBed:
            raise ValueError("File needs to be in bigBed format!")

    def get_filters(self):
        """
        Return filters

        Returns
        -------

        """
        return self._filters

    def set_filters(self, key, value):
        """
        Set filter, records matching filters will be labeled in the track

        Parameters
        ----------
        key : str
            Only "MAF" and "ID" are supported currently
        value : numeric
            min value
        Returns
        -------

        """
        if key in self._filter_supported_fields:
            self._filters[key] = value

    def _get(self, chromosome, start, end):
        entries = []
        try:
            for entry in self._bb.entries(chromosome, start, end):
                entries.append(entry)
        except:
            pass
        return entries

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw lollipop plot for mutations

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
        super(UCSCMutationTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        entries = self._get(chromosome=chromosome, start=start, end=end)
        xs = []
        ys = []
        backward_mapping = dict()
        try:
            for mutation_info in entries:
                items = mutation_info[2].split("\t")
                if items[6] != "":  # require SNPs to have MAFs
                    try:
                        maf = np.mean(
                            list(
                                map(
                                    float,
                                    filter(
                                        lambda x: x != "-inf" and x != "",
                                        items[6].split(","),
                                    ),
                                )
                            )
                        )
                    except Exception:
                        continue
                    # plot marker
                    ys.append(maf)
                    xs.append(mutation_info[0])
                    backward_mapping[items[0]] = (mutation_info[0], maf)
        except:
            pass

        if self.apply_color_gradient:
            norm = self.normalizer(vmin=np.min(ys), vmax=np.quantile(ys, 0.98))
            self._ax.scatter(
                xs, ys, marker="v", color=self._color_map(norm(ys)), clip_on=False
            )
        else:
            if "ID" in self._filters:
                colors = [
                    "gray",
                ] * len(xs)
                for highlight_id in self._filters["ID"]:
                    if highlight_id in backward_mapping:
                        hit_x = backward_mapping[highlight_id][0]
                        colors[xs.index(hit_x)] = "red"
                self._ax.scatter(xs, ys, marker="v", color=colors, zorder=50)
            else:
                self._ax.scatter(xs, ys, marker="v", color="red", zorder=50)
        for x, y in zip(xs, ys):
            self._ax.add_line(
                Line2D(
                    (x, x), (0, y), color=self.line_color, linewidth=self.line_width
                )
            )

        if "ID" in self._filters:
            texts = []
            for highlight_id in self._filters["ID"]:
                if highlight_id in backward_mapping:
                    texts.append(
                        self._ax.text(
                            backward_mapping[highlight_id][0],
                            backward_mapping[highlight_id][1],
                            highlight_id,
                        )
                    )
            try:
                from adjustText import adjust_text

                adjust_text(texts)
            except ImportError:
                pass


class BigBed6Track(BedTrack):
    """Standard BigBed6 track."""

    _bb: Any = PrivateAttr(default=None)
    _bb_obj: Any = PrivateAttr(default=None)
    _BigBedRecord: Any = PrivateAttr(default=None)
    _filters: dict = PrivateAttr(default_factory=dict)
    _filter_supported_fields: set = PrivateAttr(
        default_factory=lambda: {
            "contig",
            "start",
            "end",
            "name",
            "score",
            "strand",
        }
    )

    def _get(self, chromosome, start, end):
        results = []
        n_expected_fields = len(self._fields)
        warning = 0
        try:
            for entry in self._bb_obj.entries(chromosome, start, end):
                row = [chromosome, entry[0], entry[1]]
                other_items = entry[2].strip().split("\t")
                row.extend(other_items)
                if len(row) > n_expected_fields:
                    warning = len(row)
                    row = row[:n_expected_fields]
                results.append(self._BigBedRecord._make(row))
            if warning > 0:
                warnings.warn(
                    f"Input bigBed should only have {n_expected_fields} fields "
                    f"while {warning} fields are observed. Only the first {n_expected_fields} are used.",
                    RuntimeWarning,
                )
            return sorted(results, key=lambda x: x.start)
        except (ValueError, TypeError):
            return []

    def _open_source(self) -> None:
        self._bed_file = self.track
        self._fields = ("contig", "start", "end", "name", "score", "strand")
        if not os.path.exists(self.track) and not self.track.startswith("http"):
            raise ValueError
        self._bb = pyBigWig.open(self.track)
        if not self._bb.isBigBed:
            raise ValueError("File needs to be in bigBed format!")
        self._bb_obj = pyBigWig.open(self.track)
        self._BigBedRecord = namedtuple("BigBedRecord", self._fields)
        self._rgb_check = re.compile(r"(\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})")
        self._small_relative = 0
        if self.plot_thickness is None:
            self.plot_thickness = False
        if self.color is None:
            self.color = "#A1A1A1"
        if self.edge_color is None:
            self.edge_color = "#6E6E6E"

    def get_filters(self):
        """
        Return filters

        Returns
        -------

        """
        return self._filters

    def set_filters(self, key, value):
        """
        Set filter, records with matching names will be labeled

        Parameters
        ----------
        key : str
            Currently, only `name` is supported
        value : str

        Returns
        -------

        """
        if key in self._filter_supported_fields:
            self._filters[key] = value

    # def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
    #     """
    #     Draw BB6 track
    #
    #     Parameters
    #     ----------
    #     chromosome : str
    #         name of the chromosome/contig
    #     start : int
    #         start of the ROI/window, 0-based
    #     end : int
    #         end of the ROI/window, 0-based
    #     ax : :class:`matplotlib.pyplot.Axes`
    #         matplotlib.pyplot.Axes for this track
    #     index : int
    #         The first subplot (track), index==0, will have its top border and xticks shown up
    #     kwargs :
    #
    #     Returns
    #     -------
    #
    #     """
    #     super(BigBed6Track, self)._draw_track(chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs)
    #     import matplotlib.pyplot as plt
    #     fig = plt.gcf()
    #     bbox = self._ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    #     self._ax.set_xlim((start, end))
    #
    #     self.small_relative = 0.004 * (end - start)
    #     backward_mapping = defaultdict(list)
    #     for lane in self._lane_registries:
    #         for interval in lane.features:
    #             start_loc = int(interval.start)
    #             end_loc = int(interval.end)
    #             visible_start = max(start_loc, start)
    #             visible_end = min(end_loc, end)
    #             active_lane = lane.offset
    #
    #             real_active_line = (self._patch_height + self._lane_space) * active_lane
    #
    #             backward_mapping[interval.name].append((start_loc, -1 * real_active_line - (self._patch_height / 2)))
    #             is_highlight = 1 if "name" in self._filters and interval.name in self._filters[
    #                 "name"] else 0
    #             rec = Rectangle(xy=(start_loc, -1 * real_active_line - (self._patch_height / 2)),
    #                             width=end_loc - start_loc,
    #                             height=self._patch_height,
    #                             facecolor="red" if is_highlight else self.color,
    #                             alpha=1 if is_highlight else 0.5)
    #             self._ax.add_patch(rec)
    #
    #             if "name" in dir(interval) and self.show_name:
    #                 if start_loc > start and interval.strand == "+":
    #                     self._ax.text(x=start_loc - self.small_relative, y=-1 * real_active_line,
    #                                   color=self._font_color, size=self._font_size,
    #                                   s=interval.name, ha="right", va="center", clip_on=False, zorder=101)
    #                 elif end_loc < end and interval.strand == "-":
    #                     self._ax.text(x=end_loc + self.small_relative, y=-1 * real_active_line,
    #                                   color=self._font_color, size=self._font_size,
    #                                   s=interval.name, ha="left", va="center", clip_on=False, zorder=101)
    #                 else:
    #                     self._ax.text(x=(visible_end + visible_start) / 2, y=-1 * real_active_line,
    #                                   color=self._font_color, size=self._font_size,
    #                                   s=interval.name, ha="center", va="center", clip_on=False,
    #                                   bbox=dict(boxstyle="round", fc="w",
    #                                             alpha=self._font_box_alpha, lw=0.1), zorder=101
    #                                   )
    #
    #     self._ax.set_yticks([])
    #     # remove minor ticks
    #     self._ax.set_yticks([], minor=True)
    #
    #     n = len(self._lane_registries)
    #     self._ax.set_ylim((-1 * n, 1.5))
    #
    #     if "name" in self._filters:
    #         texts = []
    #         for highlight_name in self._filters["name"]:
    #             if highlight_name in backward_mapping:
    #                 for hl in backward_mapping[highlight_name]:
    #                     texts.append(self._ax.text(hl[0],
    #                                                hl[1],
    #                                                highlight_name))
    #         try:
    #             from adjustText import adjust_text
    #             adjust_text(texts, arrowprops=dict(arrowstyle='-', color=self.font_color))
    #         except ImportError:
    #             pass
    #
    #     if index != 0:
    #         # remove major ticks
    #         self._ax.set_xticks([])
    #         # remove minor ticks
    #         self._ax.set_xticks([], minor=True)
    #         # self.ax.margins(0)
