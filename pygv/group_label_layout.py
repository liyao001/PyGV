"""Automatic placement of group labels to the left of track decorations."""

from contextlib import nullcontext
from warnings import warn

from matplotlib.layout_engine import LayoutEngine, PlaceHolderLayoutEngine

# Physical gap between group text, connecting lines, and track decorations.
_PAD_POINTS = 4.0
_LINE_WIDTH_POINTS = 1.0
_MAX_LAYOUT_ITERS = 6
_MIN_SIGNAL_FRAC = 0.30
_MIN_SIGNAL_INCHES = 0.75
_LEFT_TOL = 1e-4
_OVERLAP_EPS_PX = 1.0


def snapshot_subplotpars(fig):
    """Return the figure's current subplot parameters as a dict."""
    sp = fig.subplotpars
    return {
        "left": sp.left,
        "right": sp.right,
        "top": sp.top,
        "bottom": sp.bottom,
        "wspace": sp.wspace,
        "hspace": sp.hspace,
    }


def wrap_layout_engine(existing):
    """Return an inner engine that should run before the group-label gutter."""
    if existing is None or isinstance(existing, PlaceHolderLayoutEngine):
        return None
    if isinstance(existing, GroupLabelLayoutEngine):
        return existing._inner
    return existing


class GroupLabelLayoutEngine(LayoutEngine):
    """Pack group labels into a left gutter after the underlying figure layout."""

    _adjust_compatible = True
    _colorbar_gridspec = True

    def __init__(self, viewer, axes, artists, inner=None, base_subplots=None):
        super().__init__()
        if inner is not None:
            self._adjust_compatible = inner.adjust_compatible
            self._colorbar_gridspec = inner.colorbar_gridspec
        self._viewer = viewer
        self._axes = list(axes)
        self._artists = artists
        self._inner = inner
        self._base_subplots = dict(base_subplots or {})
        self._warned_unfit = False

    def set(self, **kwargs):
        self._params.update(kwargs)

    def execute(self, fig):
        artists = []
        for item in self._artists:
            artists.append(item["text"])
            artists.append(item["line"])

        if self._inner is not None:
            vis = [artist.get_visible() for artist in artists]
            for artist in artists:
                artist.set_visible(False)
            try:
                self._inner.execute(fig)
            except ValueError:
                pass
            for artist, visible in zip(artists, vis):
                artist.set_visible(visible)
        elif "left" in self._base_subplots:
            # Restore only left so later top/right adjustments (e.g. axis
            # marks) survive redraws, while the group gutter does not grow.
            fig.subplots_adjust(left=self._base_subplots["left"])

        _layout_group_labels(fig, self._axes, self._artists, warned=self)


def _layout_group_labels(fig, axes, artists, warned=None):
    renderer = fig._get_renderer()
    draw_disabled = getattr(renderer, "_draw_disabled", nullcontext)
    with draw_disabled():
        auto_items = [item for item in artists if item["config"].x is None]
        if auto_items:
            for _ in range(_MAX_LAYOUT_ITERS):
                _prepare_axis_decorations(axes, renderer)
                left_before = fig.subplotpars.left
                _adjust_left_for_automatic_labels(
                    fig, axes, auto_items, renderer, warned=warned
                )
                if abs(fig.subplotpars.left - left_before) < _LEFT_TOL:
                    break
            _prepare_axis_decorations(axes, renderer)
            _place_automatic_labels(fig, axes, auto_items, renderer)
        _update_manual_and_y(fig, axes, artists, renderer)


def _prepare_axis_decorations(axes, renderer):
    """Update tick, label, and offset-text positions before measuring."""
    for ax in axes:
        if not ax.get_visible():
            continue
        ax.yaxis._update_ticks()
        ax.get_tightbbox(renderer, call_axes_locator=False, for_layout_only=True)


def _update_manual_and_y(fig, axes, artists, renderer=None):
    """Apply manual x formulas and keep every line aligned to its tracks."""
    n_axes = len(axes)
    for item in artists:
        config = item["config"]
        if not _config_axes_in_range(config, n_axes):
            continue
        y0, y1 = _group_y_span_figure(axes, config)
        y_center = (y0 + y1) / 2.0
        if config.x is not None:
            x = float(config.x)
            x_line_offset = float(config.x_line_offset)
            item["text"].set_x(x - x_line_offset)
            item["text"].set_y(y_center)
            item["text"].set_ha("center")
            item["line"].set_data([x + x_line_offset, x + x_line_offset], [y0, y1])
        else:
            if renderer is not None:
                _, text_height = _text_display_size(item["text"], renderer)
                y_center = _clamp_text_y(fig, y_center, text_height, renderer)
            item["text"].set_y(y_center)
            xdata = item["line"].get_xdata()
            if len(xdata) >= 1:
                item["line"].set_data([xdata[0], xdata[0]], [y0, y1])


def _config_axes_in_range(config, n_axes):
    return (
        0 <= config.start_track_idx < n_axes
        and 0 <= config.end_track_idx < n_axes
    )


def _group_y_span_figure(axes, config):
    bbox_top = axes[config.start_track_idx].get_position()
    bbox_bottom = axes[config.end_track_idx].get_position()
    return bbox_bottom.y0, bbox_top.y1


def _column_plan(fig, axes, auto_items, renderer):
    pad_px = renderer.points_to_pixels(_PAD_POINTS)
    line_px = max(renderer.points_to_pixels(_LINE_WIDTH_POINTS), 1.0)
    decoration_left = _leftmost_decoration_display_x(axes, renderer)
    columns = _assign_columns(fig, axes, auto_items, renderer)

    cursor = decoration_left
    placements = {}
    for column in columns:
        max_text_width = max(item["text_width"] for item in column)
        cursor -= pad_px
        line_x = cursor - line_px / 2.0
        cursor -= line_px
        cursor -= pad_px
        for item in column:
            text_x = cursor - item["text_width"] / 2.0
            placements[id(item["text"])] = (text_x, line_x, item["text_width"])
        cursor -= max_text_width
    return placements, cursor, pad_px


def _adjust_left_for_automatic_labels(fig, axes, auto_items, renderer, warned=None):
    _, leftmost, pad_px = _column_plan(fig, axes, auto_items, renderer)
    fig_width_px = fig.bbox.width
    if fig_width_px <= 0:
        return
    needed_left_shift_px = pad_px - leftmost if leftmost < pad_px else 0.0
    if needed_left_shift_px <= 0:
        return
    needed_left = fig.subplotpars.left + needed_left_shift_px / fig_width_px
    max_left = _max_left_for_signal(fig)
    if needed_left > max_left + _LEFT_TOL:
        needed_left = max_left
        _warn_unfit(warned)
    if needed_left > fig.subplotpars.left + _LEFT_TOL:
        fig.subplots_adjust(left=needed_left)


def _place_automatic_labels(fig, axes, auto_items, renderer):
    placements, _, _ = _column_plan(fig, axes, auto_items, renderer)
    inv = fig.transFigure.inverted()
    for item in auto_items:
        text_x, line_x, text_width = placements[id(item["text"])]
        text_x, line_x = _clamp_x_inside_figure(
            fig, text_x, line_x, text_width, renderer
        )
        fx, _ = inv.transform((text_x, 0.0))
        lx, _ = inv.transform((line_x, 0.0))
        y0, y1 = _group_y_span_figure(axes, item["config"])
        _, text_height = _text_display_size(item["text"], renderer)
        y_center = _clamp_text_y(fig, (y0 + y1) / 2.0, text_height, renderer)
        item["text"].set_x(fx)
        item["text"].set_y(y_center)
        item["text"].set_ha("center")
        item["line"].set_data([lx, lx], [y0, y1])


def _clamp_text_y(fig, y_center, text_height, renderer):
    """Shift vertical text just enough to stay inside the figure."""
    pad_px = renderer.points_to_pixels(_PAD_POINTS)
    fig_height_px = fig.bbox.height
    _, y_disp = fig.transFigure.transform((0.0, y_center))
    half = text_height / 2.0
    lo = pad_px + half
    hi = fig_height_px - pad_px - half
    if lo > hi:
        y_disp = fig_height_px / 2.0
    else:
        y_disp = min(max(y_disp, lo), hi)
    _, y_fig = fig.transFigure.inverted().transform((0.0, y_disp))
    return y_fig


def _clamp_x_inside_figure(fig, text_x, line_x, text_width, renderer):
    """Keep automatic text inside the figure; may overlap decorations if unfit."""
    pad_px = renderer.points_to_pixels(_PAD_POINTS)
    fig_width_px = fig.bbox.width
    half = text_width / 2.0
    lo = pad_px + half
    hi = fig_width_px - pad_px - half
    if lo > hi:
        new_text_x = fig_width_px / 2.0
    else:
        new_text_x = min(max(text_x, lo), hi)
    shift = new_text_x - text_x
    return new_text_x, line_x + shift


def _max_left_for_signal(fig):
    fig_width = fig.get_figwidth()
    min_width = max(_MIN_SIGNAL_FRAC, _MIN_SIGNAL_INCHES / max(fig_width, 1e-6))
    min_width = min(min_width, 0.60)
    return max(fig.subplotpars.right - min_width, fig.subplotpars.left)


def _warn_unfit(warned):
    if warned is not None and getattr(warned, "_warned_unfit", False):
        return
    warn(
        "The figure is too narrow to place group labels beside the track "
        "decorations while keeping a usable signal area. Use a larger "
        "fig_width or smaller labels.",
        RuntimeWarning,
        stacklevel=2,
    )
    if warned is not None:
        warned._warned_unfit = True


def _assign_columns(fig, axes, auto_items, renderer):
    measured = []
    for item in auto_items:
        width, height = _text_display_size(item["text"], renderer)
        interval = _group_interval_display(
            fig, axes, item["config"], height, renderer
        )
        measured.append(
            {
                "config": item["config"],
                "text": item["text"],
                "line": item["line"],
                "text_width": width,
                "text_height": height,
                "interval": interval,
            }
        )
    measured.sort(key=lambda m: m["config"].start_track_idx)

    columns = []
    for item in measured:
        placed = False
        for column in columns:
            if all(
                not _intervals_overlap(item["interval"], other["interval"])
                for other in column
            ):
                column.append(item)
                placed = True
                break
        if not placed:
            columns.append([item])
    return columns


def _intervals_overlap(a, b):
    lo_a, hi_a = a
    lo_b, hi_b = b
    return hi_a > lo_b + _OVERLAP_EPS_PX and hi_b > lo_a + _OVERLAP_EPS_PX


def _group_interval_display(fig, axes, config, text_height, renderer):
    y0, y1 = _group_y_span_figure(axes, config)
    _, yd0 = fig.transFigure.transform((0.0, y0))
    _, yd1 = fig.transFigure.transform((0.0, y1))
    line_lo, line_hi = (yd0, yd1) if yd0 <= yd1 else (yd1, yd0)
    y_center_fig = 0.5 * (y0 + y1)
    y_center_fig = _clamp_text_y(fig, y_center_fig, text_height, renderer)
    _, y_center = fig.transFigure.transform((0.0, y_center_fig))
    text_lo = y_center - text_height / 2.0
    text_hi = y_center + text_height / 2.0
    return min(line_lo, text_lo), max(line_hi, text_hi)


def _text_display_size(text, renderer):
    bbox = text.get_window_extent(renderer)
    return bbox.width, bbox.height


def _leftmost_decoration_display_x(axes, renderer):
    xs = []
    for ax in axes:
        if not ax.get_visible():
            continue
        ax_bbox = ax.get_window_extent(renderer)
        xs.append(ax_bbox.x0)
        for artist in _iter_yaxis_decoration_texts(ax):
            bbox = artist.get_window_extent(renderer)
            xs.append(bbox.x0)
    if not xs:
        return 0.0
    return min(xs)


def _yaxis_left_tick_labels_on(ax, which):
    """True if left-side tick labels for *which* ('major'/'minor') are enabled.

    Avoids ``Axis.get_tick_params`` (Matplotlib 3.7+) so Matplotlib 3.6 works.
    """
    if which == "minor":
        ticks = ax.yaxis.get_minor_ticks()
    else:
        ticks = ax.yaxis.get_major_ticks()
    if not ticks:
        return False
    tick = ticks[0]
    if hasattr(tick, "label1On"):
        return bool(tick.label1On)
    label = getattr(tick, "label1", None)
    if label is None:
        return False
    return bool(label.get_visible())


def _iter_yaxis_decoration_texts(ax):
    if not ax.yaxis.get_visible():
        return

    ylabel = ax.yaxis.label
    if ylabel.get_visible() and str(ylabel.get_text()).strip():
        yield ylabel

    for which in ("major", "minor"):
        if not _yaxis_left_tick_labels_on(ax, which):
            continue
        for label in ax.yaxis.get_ticklabels(which=which):
            if label.get_visible() and str(label.get_text()).strip():
                yield label

    offset = ax.yaxis.get_offset_text()
    if offset.get_visible() and str(offset.get_text()).strip():
        yield offset
