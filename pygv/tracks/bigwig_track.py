import numpy as np
import pyBigWig
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from pygv.utils import check_accessibility

from .track import NumericalTrack


class BigWigTrack(NumericalTrack):
    """
    Generic BigWig track

    Examples
    --------

    .. plot:: ../examples/plot_bigwig.py
    """

    def __init__(self, track: str, plot_type: str = "line", **kwargs):
        """

        Parameters
        ----------
        track : str or list[str]
            Path to the bigwig file(s), can be local path(s) or url(s) pointing to a remote file
        kwargs

        """
        super(BigWigTrack, self).__init__(**kwargs)
        if isinstance(track, str):
            check_accessibility(track, allow_remote=True)
            self.bw = [
                pyBigWig.open(track),
            ]
        else:
            self.bw = []
            for sub_t in track:
                check_accessibility(sub_t, allow_remote=True)
                self.bw.append(pyBigWig.open(sub_t))

        self.plot_type = plot_type

    def _get(self, chromosome, start, end, nan_as_zero=True):
        """
        Get data from bigwig file among a certain window

        Parameters
        ----------
        chromosome : str
            name of the chromosome/contig
        start : int
            start of the ROI/window, 0-based
        end : int
            end of the ROI/window, 0-based
        nan_as_zero : bool
            Set NaNs as zeros

        Returns
        -------
        xvalues : :class:`numpy.ndarray`
            x values
        values : :class:`numpy.ndarray`
            y values
        """
        values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self.bw
            ]
        ).mean(axis=0)
        values = self.data_transform(values)
        if self._scale != 1:
            values *= self._scale
        xvalues = np.arange(start, end, step=1)

        if self.stat_method is not None:
            from scipy.stats import binned_statistic

            y_new, x_new, _ = binned_statistic(
                xvalues, values, statistic=self._stat_method, bins=self._n_bins
            )
            xvalues = x_new[:-1]
            values = y_new
        # find indices where consecutive values are not equal
        keep_idx = self._merge_redundant_values(xvalues, values)
        return xvalues[keep_idx], values[keep_idx]

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw track

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
        super(BigWigTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )

        if self.plot_type == "line":
            x, y = self._get(chromosome=chromosome, start=start, end=end)
            self._ax.plot(
                x, y, color=self.color, linewidth=self.line_width, alpha=self._alpha
            )
            self._ax.fill_between(
                x, y, 0, facecolor=self.color, alpha=self._alpha, lw=self.line_width
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
    """
    Overlay BigWig tracks (signals from multiple BigWig files in the same track) in a single track

    Parameters
    ----------
    tracks : list of str or tuple[str, ...]
        List of file paths or urls. Tracks will be organized with ascending zorder.
    labels : list of str or tuple[str, ...]
        Labels for each bigwig file
    palette : str or palette instances
        Name of the palette you want to use. Matplot colormaps / seaborn palettes
    colors : list of colors or None
        List of colors you want to use. If it's `None`, then colors will be extracted from ``palette``.
    legend : bool
        Enable/disable legends
    legend_kws : dict, optional
        Keyword arguments passed to matplotlib legend.

    Examples
    --------

    .. plot:: ../examples/plot_overlaying_bigwigs.py
    """

    @property
    def labels(self):
        """
        Labels for each bigwig file
        """
        return self._labels

    @labels.setter
    def labels(self, value):
        if len(value) == len(self.bws):
            self._labels = value

    @property
    def colors(self):
        """
        Colors for each bigwig file
        """
        return self._colors

    @colors.setter
    def colors(self, value):
        if len(value) == len(self.bws):
            self._colors = value

    def __init__(
        self,
        tracks,
        labels,
        palette="Set1",
        colors=None,
        legend=True,
        legend_kws=None,
        **kwargs,
    ):
        super(OverlayingTrack, self).__init__(**kwargs)
        self.bws = []
        for track in tracks:
            check_accessibility(track, allow_remote=True)

            self.bws.append(pyBigWig.open(track))
        self._labels = labels
        from seaborn.palettes import color_palette

        self._colors = None
        if colors is None:
            self.colors = color_palette(palette=palette, n_colors=len(tracks))
        else:
            self.colors = colors
        self._legend = legend

        _default_legend_conf = {"frameon": False}
        self._legend_kws = legend_kws if isinstance(legend_kws, dict) else {}

    def _get(self, chromosome, start, end):
        value_list = []
        xvalue_list = []

        for bw_obj in self.bws:
            xvalues = np.arange(start, end, step=1)
            values = self.convert_nan_to_num(bw_obj.values(chromosome, start, end))
            values = self.data_transform(values)
            if self._scale != 1:
                values *= self._scale

            if self._stat_method is not None:
                from scipy.stats import binned_statistic

                y_new, x_new, _ = binned_statistic(
                    xvalues, values, statistic=self._stat_method, bins=self._n_bins
                )
                xvalues = x_new[:-1]
                values = y_new
            keep_idx = self._merge_redundant_values(xvalues, values)
            value_list.append(values[keep_idx])
            xvalue_list.append(xvalues[keep_idx])
        return xvalue_list, value_list

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw stacked tracks

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
        super(OverlayingTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        xs, ys = self._get(chromosome=chromosome, start=start, end=end)
        z = 0
        for x, y, c, l in zip(xs, ys, self._colors, self._labels):
            self._ax.plot(x, y, color=c, label=l, linewidth=self.line_width, zorder=z)
            self._ax.fill_between(x, y, 0, facecolor=c, alpha=self._alpha, zorder=z)
            z += 1

        if self._legend:
            self._ax.legend(**self._legend_kws)


class PairedStrandSpecificTrack(NumericalTrack):
    """
    Paired strand-specific tracks

    Parameters
    ----------
    pl_track : str or list[str]
        File path(s) or url(s) for the positive track. When multiple tracks are provided, the mean values will be used.
    mn_track : str or list[str]
        File path(s) or url(s) for the negative track. When multiple tracks are provided, the mean values will be used.
    draw_y_independently : bool
        By default, `True`, which means the output track centers at zero
        and the positive and negatives axis have identical lengths regardless to their ranges.
        Set it as `False` if you want the y-axis to reflect the dominance of signals on a strand.
    plot_type : str
        "line": line plot
        "bar": bar plot
    pos_color : color_like
        Color for positive signals, by default, #E10600
    neg_color : color_like
        Color for negative signals, by default, #0048AC
    kwargs :

    Examples
    --------

    .. plot:: ../examples/plot_paired_strand_specific_bigwigs.py
    """

    @property
    def pos_color(self):
        """
        Color for negative signals
        """
        return self._pos_color

    @pos_color.setter
    def pos_color(self, value):
        from matplotlib.colors import is_color_like

        if is_color_like(value) or value is None:
            self._pos_color = value

    @property
    def neg_color(self):
        """
        Color for negative signals
        """
        return self._neg_color

    @neg_color.setter
    def neg_color(self, value):
        from matplotlib.colors import is_color_like

        if is_color_like(value) or value is None:
            self._neg_color = value

    def __init__(
        self,
        pl_track,
        mn_track,
        draw_y_independently=True,
        plot_type: str = "line",
        **kwargs,
    ):
        super(PairedStrandSpecificTrack, self).__init__(**kwargs)
        # color
        self._pos_color = "#E10600"
        self.pos_color = kwargs.get("pos_color", "#E10600")
        self._neg_color = "#0048AC"
        self.neg_color = kwargs.get("neg_color", "#0048AC")

        for t in (pl_track, mn_track):
            if isinstance(t, str):
                check_accessibility(t, allow_remote=True)
            else:
                for sub_t in t:
                    check_accessibility(sub_t, allow_remote=True)

        if isinstance(pl_track, str):
            self.pl_bw = [
                pyBigWig.open(pl_track),
            ]
        else:
            self.pl_bw = list([pyBigWig.open(f) for f in pl_track])
        if isinstance(mn_track, str):
            self.mn_bw = [
                pyBigWig.open(mn_track),
            ]
        else:
            self.mn_bw = list([pyBigWig.open(f) for f in mn_track])
        assert len(self.pl_bw) == len(self.mn_bw)
        self._equal_space_for_pos_neg_ranges = draw_y_independently
        self.plot_type = plot_type

    def _get(self, chromosome, start, end, nan_as_zero=True):
        pl_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self.pl_bw
            ]
        ).mean(axis=0)
        mn_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self.mn_bw
            ]
        ).mean(axis=0)
        pl_values = self.data_transform(pl_values)
        if self._scale != 1:
            pl_values *= self._scale
        mn_values[mn_values > 0] *= -1
        mn_values = self.data_transform(mn_values)
        if self._scale != 1:
            mn_values *= self._scale
        xvalues = np.arange(start, end, step=1)

        if self._stat_method is not None:
            from scipy.stats import binned_statistic

            pl_new, x_new, _ = binned_statistic(
                xvalues, pl_values, statistic=self._stat_method, bins=self._n_bins
            )
            mn_new, x_new, _ = binned_statistic(
                xvalues, mn_values, statistic=self._stat_method, bins=self._n_bins
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
        """
        Draw strand-specific tracks

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
        super(PairedStrandSpecificTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        self.is_real_number_track = 1
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
    """
    Paired strandless tracks

    Parameters
    ----------
    pl_track : str or list[str]
        File path(s) or url(s) for the positive track. When multiple tracks are provided, the mean values will be used.
    mn_track : str or list[str]
        File path(s) or url(s) for the negative track. When multiple tracks are provided, the mean values will be used.
    plot_type : str
        "line": line plot
        "bar": bar plot
    kwargs :

    Examples
    --------

    .. plot:: ../examples/plot_paired_strandless_bigwigs.py
    """

    def __init__(self, pl_track, mn_track, plot_type: str = "line", **kwargs):
        super(PairedStrandlessTrack, self).__init__(
            track=pl_track, plot_type=plot_type, **kwargs
        )
        self.pl_bw = []
        self.mn_bw = []
        for t in (pl_track, mn_track):
            if isinstance(t, str):
                check_accessibility(t, allow_remote=True)
            else:
                for sub_t in t:
                    check_accessibility(sub_t, allow_remote=True)

        if isinstance(pl_track, str):
            self.pl_bw = [
                pyBigWig.open(pl_track),
            ]
        else:
            self.pl_bw = list([pyBigWig.open(f) for f in pl_track])
        if isinstance(mn_track, str):
            self.mn_bw = [
                pyBigWig.open(mn_track),
            ]
        else:
            self.mn_bw = list([pyBigWig.open(f) for f in mn_track])
        assert len(self.pl_bw) == len(self.mn_bw)

        self.plot_type = plot_type

    def _get(self, chromosome, start, end, nan_as_zero=True):
        pl_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self.pl_bw
            ]
        ).mean(axis=0)
        mn_values = np.stack(
            [
                (
                    self.convert_nan_to_num(_bw.values(chromosome, start, end))
                    if nan_as_zero
                    else _bw.values(chromosome, start, end, numpy=True)
                )
                for _bw in self.mn_bw
            ]
        ).mean(axis=0)
        pl_values = self.data_transform(pl_values)
        if self._scale != 1:
            pl_values *= self._scale
        mn_values = np.abs(mn_values)
        mn_values = self.data_transform(mn_values)
        if self._scale != 1:
            mn_values *= self._scale
        xvalues = np.arange(start, end, step=1)
        y_values = np.nansum(np.stack([pl_values, mn_values]), axis=0)

        if self._stat_method is not None:
            from scipy.stats import binned_statistic

            y_new, x_new, _ = binned_statistic(
                xvalues, y_values, statistic=self._stat_method, bins=self._n_bins
            )
            xvalues = x_new[:-1]
            y_values = y_new
        keep_idx = self._merge_redundant_values(xvalues, y_values)
        return xvalues[keep_idx], y_values[keep_idx]
