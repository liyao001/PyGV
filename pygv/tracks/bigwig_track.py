from __future__ import annotations

from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import pyBigWig
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from pydantic import AliasChoices, Field, PrivateAttr, model_validator
from seaborn.palettes import color_palette

from pygv.utils import check_accessibility

from .track import NumericalTrack
from .types import Color, PlotType


class BigWigTrack(NumericalTrack):
    """Generic BigWig track."""

    track: Union[str, List[str]] = Field(
        description="Path to the bigwig file(s), local path(s) or URL(s)"
    )
    plot_type: PlotType = Field(
        default="line",
        kw_only=True,
        description='Plot style: "line" or "bar"',
    )

    _bw: list = PrivateAttr(default_factory=list)

    def __init__(self, track: Union[str, List[str]], **data: Any) -> None:
        super().__init__(track=track, **data)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        self._open_bigwigs()

    def _open_bigwigs(self) -> None:
        if isinstance(self.track, str):
            check_accessibility(self.track, allow_remote=True)
            self._bw = [pyBigWig.open(self.track)]
        else:
            opened = []
            for sub_t in self.track:
                check_accessibility(sub_t, allow_remote=True)
                opened.append(pyBigWig.open(sub_t))
            self._bw = opened

    def _get(self, chromosome, start, end, nan_as_zero=True):
        values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self._bw
            ]
        ).mean(axis=0)
        values = self.data_transform(values)
        if self.scale != 1:
            values *= self.scale
        xvalues = np.arange(start, end, step=1)

        if self.stat_method is not None:
            from scipy.stats import binned_statistic

            y_new, x_new, _ = binned_statistic(
                xvalues, values, statistic=self.stat_method, bins=self.n_bins
            )
            xvalues = x_new[:-1]
            values = y_new
        keep_idx = self._merge_redundant_values(xvalues, values)
        return xvalues[keep_idx], values[keep_idx]

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )

        if self.plot_type == "line":
            x, y = self._get(chromosome=chromosome, start=start, end=end)
            self._ax.plot(
                x, y, color=self.color, linewidth=self.line_width, alpha=self.alpha
            )
            self._ax.fill_between(
                x, y, 0, facecolor=self.color, alpha=self.alpha, lw=self.line_width
            )
        elif self.plot_type == "bar":
            X, Y = self._get(
                chromosome=chromosome, start=start, end=end, nan_as_zero=False
            )
            w_rect = 1 / (end - start)
            half_w = w_rect / 2.0
            rects = []
            for x, y in zip(X, Y):
                if not np.isnan(y):
                    rec_y_start = 0 if y > 0 else y
                    rec_y_end = np.abs(y)
                    rec = Rectangle(
                        xy=(x - half_w, rec_y_start),
                        width=w_rect,
                        height=rec_y_end,
                        edgecolor=self.edge_color,
                        facecolor=self.color,
                        linewidth=self.line_width,
                    )
                    rects.append(rec)
            self._ax.scatter(X, Y, s=0)
            self._ax.add_collection(
                PatchCollection(
                    rects,
                    edgecolors=self.edge_color,
                    facecolors=self.color,
                    linewidths=self.line_width,
                    zorder=100,
                    clip_on=True,
                )
            )


class OverlayingTrack(NumericalTrack):
    """Overlay BigWig tracks (signals from multiple files) in a single track."""

    tracks: Union[List[str], Tuple[str, ...]] = Field(
        description="File paths or urls. Tracks are organized with ascending zorder."
    )
    labels: Union[List[str], Tuple[str, ...]] = Field(
        description="Labels for each bigwig file"
    )
    palette: str = Field(
        default="Set1",
        kw_only=True,
        description="Name of the palette (matplotlib / seaborn)",
    )
    colors: Optional[Sequence] = Field(
        default=None,
        kw_only=True,
        description="Colors for each file. If None, colors are taken from palette.",
    )
    legend: bool = Field(default=True, kw_only=True, description="Enable/disable legends")
    legend_kws: Optional[dict] = Field(
        default=None,
        kw_only=True,
        description="Keyword arguments passed to matplotlib legend",
    )

    _bws: list = PrivateAttr(default_factory=list)
    _legend_kws: dict = PrivateAttr(default_factory=dict)

    def __init__(
        self,
        tracks: Union[List[str], Tuple[str, ...]],
        labels: Union[List[str], Tuple[str, ...]],
        **data: Any,
    ) -> None:
        super().__init__(tracks=tracks, labels=labels, **data)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        opened = []
        for track in self.tracks:
            check_accessibility(track, allow_remote=True)
            opened.append(pyBigWig.open(track))
        self._bws = opened
        if len(self.labels) != len(self._bws):
            raise ValueError("labels must have the same length as tracks")
        if self.colors is None:
            self.colors = color_palette(palette=self.palette, n_colors=len(self.tracks))
        elif len(self.colors) != len(self._bws):
            raise ValueError("colors must have the same length as tracks")
        self._legend_kws = self.legend_kws if isinstance(self.legend_kws, dict) else {}

    def _get(self, chromosome, start, end):
        value_list = []
        xvalue_list = []

        for bw_obj in self._bws:
            xvalues = np.arange(start, end, step=1)
            values = self.convert_nan_to_num(bw_obj.values(chromosome, start, end))
            values = self.data_transform(values)
            if self.scale != 1:
                values *= self.scale

            if self.stat_method is not None:
                from scipy.stats import binned_statistic

                y_new, x_new, _ = binned_statistic(
                    xvalues, values, statistic=self.stat_method, bins=self.n_bins
                )
                xvalues = x_new[:-1]
                values = y_new
            keep_idx = self._merge_redundant_values(xvalues, values)
            value_list.append(values[keep_idx])
            xvalue_list.append(xvalues[keep_idx])
        return xvalue_list, value_list

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        xs, ys = self._get(chromosome=chromosome, start=start, end=end)
        z = 0
        for x, y, c, l in zip(xs, ys, self.colors, self.labels):
            self._ax.plot(x, y, color=c, label=l, linewidth=self.line_width, zorder=z)
            self._ax.fill_between(x, y, 0, facecolor=c, alpha=self.alpha, zorder=z)
            z += 1

        if self.legend:
            self._ax.legend(**self._legend_kws)


class PairedStrandSpecificTrack(NumericalTrack):
    """Paired strand-specific BigWig tracks."""

    pl_track: Union[str, List[str]] = Field(
        description="File path(s) or url(s) for the positive track"
    )
    mn_track: Union[str, List[str]] = Field(
        description="File path(s) or url(s) for the negative track"
    )
    plot_type: PlotType = Field(default="line", kw_only=True)
    pos_color: Color = Field(
        default="#E10600", kw_only=True, description="Color for positive signals"
    )
    neg_color: Color = Field(
        default="#0048AC", kw_only=True, description="Color for negative signals"
    )
    equal_space_for_pos_neg_ranges: bool = Field(
        default=True,
        kw_only=True,
        validation_alias=AliasChoices(
            "equal_space_for_pos_neg_ranges", "draw_y_independently"
        ),
        description="Center at zero with equal positive/negative axis lengths",
    )

    _pl_bw: list = PrivateAttr(default_factory=list)
    _mn_bw: list = PrivateAttr(default_factory=list)

    def __init__(
        self,
        pl_track: Union[str, List[str]],
        mn_track: Union[str, List[str]],
        **data: Any,
    ) -> None:
        super().__init__(pl_track=pl_track, mn_track=mn_track, **data)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        for t in (self.pl_track, self.mn_track):
            if isinstance(t, str):
                check_accessibility(t, allow_remote=True)
            else:
                for sub_t in t:
                    check_accessibility(sub_t, allow_remote=True)

        if isinstance(self.pl_track, str):
            self._pl_bw = [pyBigWig.open(self.pl_track)]
        else:
            self._pl_bw = [pyBigWig.open(f) for f in self.pl_track]
        if isinstance(self.mn_track, str):
            self._mn_bw = [pyBigWig.open(self.mn_track)]
        else:
            self._mn_bw = [pyBigWig.open(f) for f in self.mn_track]
        if len(self._pl_bw) != len(self._mn_bw):
            raise ValueError("pl_track and mn_track must have the same number of files")

    def _get(self, chromosome, start, end, nan_as_zero=True):
        pl_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self._pl_bw
            ]
        ).mean(axis=0)
        mn_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self._mn_bw
            ]
        ).mean(axis=0)
        pl_values = self.data_transform(pl_values)
        if self.scale != 1:
            pl_values *= self.scale
        mn_values[mn_values > 0] *= -1
        mn_values = self.data_transform(mn_values)
        if self.scale != 1:
            mn_values *= self.scale
        xvalues = np.arange(start, end, step=1)

        if self.stat_method is not None:
            from scipy.stats import binned_statistic

            pl_new, x_new, _ = binned_statistic(
                xvalues, pl_values, statistic=self.stat_method, bins=self.n_bins
            )
            mn_new, x_new, _ = binned_statistic(
                xvalues, mn_values, statistic=self.stat_method, bins=self.n_bins
            )
            xvalues = x_new[:-1]
            pl_values = pl_new
            mn_values = mn_new
        pl_keep_idx = self._merge_redundant_values(xvalues, pl_values)
        mn_keep_idx = self._merge_redundant_values(xvalues, mn_values)
        return (
            xvalues[pl_keep_idx],
            xvalues[mn_keep_idx],
            pl_values[pl_keep_idx],
            mn_values[mn_keep_idx],
        )

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        self._is_real_number_track = 1
        if self.plot_type == "line":
            x_plus, x_minus, y_plus, y_minus = self._get(
                chromosome=chromosome, start=start, end=end
            )
            self._ax.axhline(0, color="black", linewidth=0.5)
            self._ax.plot(
                x_plus,
                y_plus,
                color=self.pos_color,
                alpha=self.alpha,
                linewidth=self.line_width,
            )
            self._ax.fill_between(
                x_plus, y_plus, 0, facecolor=self.pos_color, alpha=self.alpha
            )
            self._ax.plot(
                x_minus,
                y_minus,
                color=self.neg_color,
                alpha=self.alpha,
                linewidth=self.line_width,
            )
            self._ax.fill_between(
                x_minus, y_minus, 0, facecolor=self.neg_color, alpha=self.alpha
            )
        elif self.plot_type == "bar":
            X_plus, X_minus, Y_plus, Y_minus = self._get(
                chromosome=chromosome, start=start, end=end, nan_as_zero=False
            )
            w_rect = 1 / (end - start)
            half_w = w_rect / 2.0

            for X, Y, color in (
                (X_plus, Y_plus, self.pos_color),
                (X_minus, Y_minus, self.neg_color),
            ):
                rects = []
                for x, y in zip(X, Y):
                    if not np.isnan(y):
                        rec_y_start = 0 if y > 0 else y
                        rec_y_end = np.abs(y)
                        rec = Rectangle(
                            xy=(x - half_w, rec_y_start),
                            width=w_rect,
                            height=rec_y_end,
                            edgecolor=color,
                            facecolor=color,
                            linewidth=self.line_width,
                        )
                        rects.append(rec)
                self._ax.add_collection(
                    PatchCollection(
                        rects,
                        edgecolors=color,
                        facecolors=color,
                        linewidths=self.line_width,
                        zorder=100,
                        clip_on=True,
                    )
                )
            self._ax.scatter(X_plus, Y_plus, s=0)
            self._ax.scatter(X_minus, Y_minus, s=0)


class PairedStrandSpecificTracks(PairedStrandSpecificTrack):
    pass


class PairedStrandlessTrack(BigWigTrack):
    """Paired strandless tracks (plus and minus signals summed as positive)."""

    track: Union[str, List[str]] = Field(default="", kw_only=True)
    pl_track: Union[str, List[str]] = Field(
        description="File path(s) or url(s) for the positive track"
    )
    mn_track: Union[str, List[str]] = Field(
        description="File path(s) or url(s) for the negative track"
    )

    _pl_bw: list = PrivateAttr(default_factory=list)
    _mn_bw: list = PrivateAttr(default_factory=list)

    def __init__(
        self,
        pl_track: Union[str, List[str]],
        mn_track: Union[str, List[str]],
        **data: Any,
    ) -> None:
        super().__init__(
            track=pl_track, pl_track=pl_track, mn_track=mn_track, **data
        )

    @model_validator(mode="before")
    @classmethod
    def _fill_track_from_pl(cls, data):
        if isinstance(data, dict) and not data.get("track"):
            data = dict(data)
            data["track"] = data.get("pl_track") or ""
        return data

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        for t in (self.pl_track, self.mn_track):
            if isinstance(t, str):
                check_accessibility(t, allow_remote=True)
            else:
                for sub_t in t:
                    check_accessibility(sub_t, allow_remote=True)

        if isinstance(self.pl_track, str):
            self._pl_bw = [pyBigWig.open(self.pl_track)]
        else:
            self._pl_bw = [pyBigWig.open(f) for f in self.pl_track]
        if isinstance(self.mn_track, str):
            self._mn_bw = [pyBigWig.open(self.mn_track)]
        else:
            self._mn_bw = [pyBigWig.open(f) for f in self.mn_track]
        if len(self._pl_bw) != len(self._mn_bw):
            raise ValueError("pl_track and mn_track must have the same number of files")

    def _get(self, chromosome, start, end, nan_as_zero=True):
        pl_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self._pl_bw
            ]
        ).mean(axis=0)
        mn_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self._mn_bw
            ]
        ).mean(axis=0)
        pl_values = self.data_transform(pl_values)
        if self.scale != 1:
            pl_values *= self.scale
        mn_values = np.abs(mn_values)
        mn_values = self.data_transform(mn_values)
        if self.scale != 1:
            mn_values *= self.scale
        xvalues = np.arange(start, end, step=1)
        y_values = np.nansum(np.stack([pl_values, mn_values]), axis=0)

        if self.stat_method is not None:
            from scipy.stats import binned_statistic

            y_new, x_new, _ = binned_statistic(
                xvalues, y_values, statistic=self.stat_method, bins=self.n_bins
            )
            xvalues = x_new[:-1]
            y_values = y_new
        keep_idx = self._merge_redundant_values(xvalues, y_values)
        return xvalues[keep_idx], y_values[keep_idx]
