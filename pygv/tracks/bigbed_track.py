import os
import warnings
from collections import namedtuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig
from matplotlib.lines import Line2D

from .bed_track import BedTrack
from .track import NumericalTrack


class UCSCMutationTrack(NumericalTrack):
    """
    Lollipop plot from UCSC-style mutational bigwig files

    Parameters
    ----------
    track :
    kwargs :
        line_color : color_like
            :attr:`line_color`
        apply_color_gradient : bool
            :attr:`apply_color_gradient`
        color_map : str or `matplotlib.pyplot.cm`
            :attr:`color_map`
    """

    @property
    def line_color(self):
        """
        Line color, by default, red
        """
        return self._line_color

    @line_color.setter
    def line_color(self, value):
        try:
            self._line_color = value
        except:
            pass

    @property
    def apply_color_gradient(self):
        """
        Set it to True to apply color gradient to marks, by default, False
        """
        return self._apply_color_gradient

    @apply_color_gradient.setter
    def apply_color_gradient(self, value):
        try:
            self._apply_color_gradient = bool(value)
        except:
            pass

    @property
    def color_map(self):
        """
        Color map for markers
        """
        return self._color_map

    @color_map.setter
    def color_map(self, value):
        try:
            self._color_map = matplotlib.colors.ListedColormap(
                value(np.linspace(0, 1, 20))[:-5, :-1]
            )
        except:
            pass

    def __init__(self, track, **kwargs):
        super(UCSCMutationTrack, self).__init__(**kwargs)
        self._filters = dict()
        self._filter_supported_fields = {"MAF", "ID"}
        self._line_color = "red"
        self.line_color = kwargs.pop("line_color", "red")
        self._apply_color_gradient = False
        self.apply_color_gradient = kwargs.get("apply_color_gradient", False)
        self._color_map = plt.cm.Reds
        self.color_map = kwargs.get("color_map", plt.cm.Reds)
        self._normalizer = kwargs.get("normalizer", matplotlib.colors.Normalize)
        if not os.path.exists(track) and not track.startswith("http"):
            raise ValueError

        self.bb = pyBigWig.open(track)
        if not self.bb.isBigBed:
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
            for entry in self.bb.entries(chromosome, start, end):
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
                    except Exception as e:
                        continue
                    # plot marker
                    ys.append(maf)
                    xs.append(mutation_info[0])
                    backward_mapping[items[0]] = (mutation_info[0], maf)
        except:
            pass

        if self._apply_color_gradient:
            norm = self._normalizer(vmin=np.min(ys), vmax=np.quantile(ys, 0.98))
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
                    (x, x), (0, y), color=self._line_color, linewidth=self.line_width
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
    """
    Standard BigBed6 track

    Parameters
    ----------
    track : str
        File path or url to a bigBed file
    kwargs : dict
        The same as :class:`pygv.tracks.track.AnnotationTrack`

    Examples
    --------

    .. plot:: ../examples/plot_bigbed.py
    """

    def _get(self, chromosome, start, end):
        results = []
        n_expected_fields = len(self.fields)
        warning = 0
        try:
            for entry in self.bb_obj.entries(chromosome, start, end):
                row = [chromosome, entry[0], entry[1]]
                other_items = entry[2].strip().split("\t")
                row.extend(other_items)
                # in cases where there are more than n_expected_fields
                # we only use the first several fields
                if len(row) > n_expected_fields:
                    warning = len(row)
                    row = row[:n_expected_fields]
                results.append(self.BigBedRecord._make(row))
            if warning > 0:
                warnings.warn(
                    f"Input bigBed should only have {n_expected_fields} fields "
                    f"while {warning} fields are observed. Only the first {n_expected_fields} are used.",
                    RuntimeWarning,
                )
            return sorted(results, key=lambda x: x.start)
        except (ValueError, TypeError):
            # in case no feature is available in that window
            return []

    def __init__(self, track, **kwargs):
        kwargs["_is_bb"] = True
        super(BigBed6Track, self).__init__(track, **kwargs)
        # parse bed file
        self.bed_file = track
        self.fields = ("contig", "start", "end", "name", "score", "strand")

        if not os.path.exists(track) and not track.startswith("http"):
            raise ValueError

        self.bb = pyBigWig.open(track)
        if not self.bb.isBigBed:
            raise ValueError("File needs to be in bigBed format!")

        self.bb_obj = pyBigWig.open(track)

        self.BigBedRecord = namedtuple("BigBedRecord", self.fields)

        self._filters = dict()
        self._filter_supported_fields = {
            "contig",
            "start",
            "end",
            "name",
            "score",
            "strand",
        }

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
