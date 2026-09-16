"""Geometry tests for automatic group-label placement."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pygv.configs.label import GroupLabelConfig
from pygv.group_label_layout import (
    _leftmost_decoration_display_x,
    _prepare_axis_decorations,
)
from pygv.tracks.track import DynamicValueTrack
from pygv.viewer import GenomeViewer

_REGION_LEN = 40
_EPS_PX = 1.5


@pytest.fixture(autouse=True)
def _cleanup_mpl():
    orig = matplotlib.rcParams.copy()
    yield
    matplotlib.rcParams.update(orig)
    plt.close("all")


def _values(scale=1.0, offset=0.0):
    return offset + scale * np.linspace(0.2, 1.0, _REGION_LEN)


def _add_dynamic_tracks(gv, names, values=None, heights=None, show_range=True, **kwargs):
    heights = heights or [0.35] * len(names)
    values = values or [_values(i + 1) for i in range(len(names))]
    tracks = []
    for i, name in enumerate(names):
        track_kwargs = dict(kwargs)
        track = DynamicValueTrack(
            "",
            name=name,
            height=heights[i],
            show_range=show_range,
            **track_kwargs,
        )
        track.values = values[i]
        gv.add_track(track)
        tracks.append(track)
    return tracks


def _plot(gv, **kwargs):
    plot_kwargs = dict(fig_width=8, height_scale_factor=0.6, hide_chromosome_name=True)
    plot_kwargs.update(kwargs)
    axs = gv.plot("chr1", 0, _REGION_LEN, **plot_kwargs)
    if not isinstance(axs, (list, tuple, np.ndarray)):
        axs = [axs]
    fig = axs[0].figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    return fig, list(axs), renderer


def _group_artists(fig):
    engine = fig.get_layout_engine()
    assert engine is not None
    return engine._artists


def _decoration_left(axs, renderer):
    _prepare_axis_decorations(axs, renderer)
    return _leftmost_decoration_display_x(axs, renderer)


def _assert_inside_figure(bbox, fig, check_y=True):
    assert bbox.x0 >= fig.bbox.x0 - _EPS_PX
    assert bbox.x1 <= fig.bbox.x1 + _EPS_PX
    if check_y:
        assert bbox.y0 >= fig.bbox.y0 - _EPS_PX
        assert bbox.y1 <= fig.bbox.y1 + _EPS_PX


def _assert_auto_geometry(fig, axs, renderer, items, require_clearance=True, check_y=True):
    deco_left = _decoration_left(axs, renderer)
    text_bboxes = []
    for item in items:
        if item["config"].x is not None:
            continue
        text_bb = item["text"].get_window_extent(renderer)
        line_bb = item["line"].get_window_extent(renderer)
        _assert_inside_figure(text_bb, fig, check_y=check_y)
        assert text_bb.x1 <= line_bb.x0 + _EPS_PX
        if require_clearance:
            assert line_bb.x1 <= deco_left + _EPS_PX
        y0, y1 = item["line"].get_data()[1]
        top = axs[item["config"].start_track_idx].get_position().y1
        bottom = axs[item["config"].end_track_idx].get_position().y0
        assert y0 == pytest.approx(bottom, abs=1e-6)
        assert y1 == pytest.approx(top, abs=1e-6)
        text_bboxes.append(text_bb)
    if require_clearance:
        for i, a in enumerate(text_bboxes):
            for b in text_bboxes[i + 1 :]:
                assert not a.overlaps(b)


def test_group_label_config_defaults_to_automatic_x():
    config = GroupLabelConfig(0, 1, "Group")
    assert config.x is None
    assert config.x_line_offset == 0.015


def test_ten_compact_tracks_middle_group_does_not_overlap_labels():
    gv = GenomeViewer(n_ticks=3, inward_ticks=True, font_size=8)
    names = [f"V{i}" for i in range(1, 11)]
    heights = [0.3] * 10
    _add_dynamic_tracks(gv, names, heights=heights)
    gv.add_group_label(3, 6, "OX2R-H3K27ac")
    fig, axs, renderer = _plot(gv)
    items = _group_artists(fig)
    assert len(items) == 1
    _assert_auto_geometry(fig, axs, renderer, items)


def test_long_names_wide_numbers_and_large_fonts():
    gv = GenomeViewer(n_ticks=3, font_size=16)
    names = [
        "upstream_coverage_replicate",
        "V2",
        "downstream_coverage_replicate",
    ]
    values = [_values(1.0), _values(12.3456), _values(2.0)]
    _add_dynamic_tracks(
        gv, names, values=values, heights=[0.5, 0.5, 0.5], y_tick_format="{:.4f}"
    )
    gv.add_group_label(0, 2, "cluster")
    fig, axs, renderer = _plot(gv, fig_width=9)
    _assert_auto_geometry(fig, axs, renderer, _group_artists(fig))


def test_hidden_range_labels_do_not_consume_gutter_space():
    def _line_x(show_range):
        gv = GenomeViewer(n_ticks=3)
        _add_dynamic_tracks(
            gv, ["", ""], values=[_values(100.0), _values(200.0)], show_range=show_range
        )
        gv.add_group_label(0, 1, "grp")
        fig, axs, renderer = _plot(gv)
        item = _group_artists(fig)[0]
        return item["line"].get_window_extent(renderer).x0, fig

    x_hidden, fig_hidden = _line_x(False)
    plt.close(fig_hidden)
    x_shown, _ = _line_x(True)
    assert x_hidden > x_shown + 1.0


def test_scientific_offset_text_is_included_in_gutter():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(
        gv,
        ["", ""],
        values=[1e6 + _values(3.0), 1e6 + _values(4.0)],
        heights=[0.5, 0.5],
    )
    gv.add_group_label(0, 1, "scaled")
    fig, axs, renderer = _plot(gv)
    for ax in axs:
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useOffset=True)
        formatter = ax.yaxis.get_major_formatter()
        if hasattr(formatter, "set_useOffset"):
            formatter.set_useOffset(True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    _prepare_axis_decorations(axs, renderer)
    offset_lefts = []
    for ax in axs:
        offset = ax.yaxis.get_offset_text()
        if offset.get_visible() and str(offset.get_text()).strip():
            offset_lefts.append(offset.get_window_extent(renderer).x0)
    assert offset_lefts, "expected scientific-notation offset text after ticklabel_format"
    deco_left = _decoration_left(axs, renderer)
    assert deco_left <= min(offset_lefts) + _EPS_PX
    _assert_auto_geometry(fig, axs, renderer, _group_artists(fig))


def test_group_autoscale_and_inward_ticks():
    gv = GenomeViewer(n_ticks=3, inward_ticks=True)
    _add_dynamic_tracks(
        gv,
        ["a", "b", "c"],
        values=[_values(1.0), _values(10.0), _values(2.0)],
        heights=[0.4, 0.8, 0.4],
    )
    gv.add_group_autoscale([0, 1])
    gv.add_group_label(0, 1, "scaled")
    fig, axs, renderer = _plot(gv)
    assert axs[0].get_ylim() == axs[1].get_ylim()
    _assert_auto_geometry(fig, axs, renderer, _group_artists(fig))


def test_overlapping_groups_use_separate_columns():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, [f"T{i}" for i in range(6)], heights=[0.4] * 6)
    gv.add_group_label(0, 3, "left-overlapping")
    gv.add_group_label(2, 5, "right-overlapping")
    fig, axs, renderer = _plot(gv, fig_width=9)
    items = _group_artists(fig)
    _assert_auto_geometry(fig, axs, renderer, items)
    line_xs = [item["line"].get_xdata()[0] for item in items]
    assert abs(line_xs[0] - line_xs[1]) > 1e-4


def test_separated_groups_share_a_column():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, [f"T{i}" for i in range(8)], heights=[0.35] * 8)
    gv.add_group_label(0, 1, "top")
    gv.add_group_label(6, 7, "bottom")
    fig, axs, renderer = _plot(gv)
    items = _group_artists(fig)
    _assert_auto_geometry(fig, axs, renderer, items)
    line_xs = [item["line"].get_xdata()[0] for item in items]
    assert line_xs[0] == pytest.approx(line_xs[1], abs=1e-6)


def test_manual_coordinates_retain_previous_formulas():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B"])
    gv.add_group_label(0, 1, "manual", x=0.02, x_line_offset=0.015)
    fig, axs, renderer = _plot(gv)
    item = _group_artists(fig)[0]
    assert item["text"].get_ha() == "center"
    assert item["text"].get_position()[0] == pytest.approx(0.005, abs=1e-6)
    assert item["line"].get_xdata()[0] == pytest.approx(0.035, abs=1e-6)
    y0, y1 = item["line"].get_data()[1]
    assert y0 == pytest.approx(axs[1].get_position().y0, abs=1e-6)
    assert y1 == pytest.approx(axs[0].get_position().y1, abs=1e-6)


def test_figures_without_groups_keep_existing_margins():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C"])
    fig, _, _ = _plot(gv)
    assert fig.get_layout_engine() is None or not hasattr(fig.get_layout_engine(), "_artists")
    left = fig.subplotpars.left

    gv2 = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv2, ["A", "B", "C"])
    fig2, _, _ = _plot(gv2)
    assert fig2.subplotpars.left == pytest.approx(left, abs=1e-6)


def test_repeated_draws_do_not_duplicate_artists_or_grow_margins():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C", "D"])
    gv.add_group_label(1, 2, "mid")
    fig, axs, renderer = _plot(gv)
    n_texts = len(fig.texts)
    n_group = len(_group_artists(fig))
    lefts = [fig.subplotpars.left]
    for _ in range(4):
        fig.canvas.draw()
        lefts.append(fig.subplotpars.left)
    assert len(fig.texts) == n_texts
    assert len(_group_artists(fig)) == n_group
    assert max(lefts) - min(lefts) < 1e-6
    _assert_auto_geometry(fig, axs, fig.canvas.get_renderer(), _group_artists(fig))


def test_resize_narrower_and_wider_recomputes_layout():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C", "D"], heights=[0.4] * 4)
    gv.add_group_label(1, 2, "mid")
    fig, axs, renderer = _plot(gv, fig_width=8)
    height = fig.get_figheight()
    fig.set_size_inches(11, height)
    fig.canvas.draw()
    _assert_auto_geometry(fig, axs, fig.canvas.get_renderer(), _group_artists(fig))
    fig.set_size_inches(6, height)
    fig.canvas.draw()
    _assert_auto_geometry(fig, axs, fig.canvas.get_renderer(), _group_artists(fig))


@pytest.mark.parametrize("force_tight_layout", [True, False])
def test_tight_layout_enabled_and_disabled(force_tight_layout):
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C"])
    gv.add_group_label(0, 2, "all")
    fig, axs, renderer = _plot(gv, force_tight_layout=force_tight_layout)
    _assert_auto_geometry(fig, axs, renderer, _group_artists(fig))


def test_png_and_pdf_export_include_group_labels(tmp_path):
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C"])
    gv.add_group_label(0, 2, "export")
    fig, axs, renderer = _plot(gv)
    items = _group_artists(fig)
    _assert_auto_geometry(fig, axs, renderer, items)
    tight = fig.get_tightbbox(renderer)
    text_inches = items[0]["text"].get_window_extent(renderer).transformed(
        fig.dpi_scale_trans.inverted()
    )
    assert text_inches.x0 >= tight.x0 - 0.05
    assert text_inches.x1 <= tight.x1 + 0.05

    png = tmp_path / "group.png"
    pdf = tmp_path / "group.pdf"
    png_tight = tmp_path / "group_tight.png"
    gv.save(str(png))
    gv.save(str(pdf))
    fig.savefig(str(png_tight), bbox_inches="tight", dpi=180)
    assert png.stat().st_size > 0
    assert pdf.stat().st_size > 0
    assert png_tight.stat().st_size > 0

    fig.set_dpi(180)
    fig.canvas.draw()
    _assert_auto_geometry(fig, axs, fig.canvas.get_renderer(), _group_artists(fig))


def test_insufficient_figure_width_warns_and_keeps_text_inside_figure():
    gv = GenomeViewer(n_ticks=3, font_size=14)
    names = [f"very_long_track_name_{i}" for i in range(6)]
    _add_dynamic_tracks(gv, names, heights=[0.45] * 6)
    gv.add_group_label(0, 3, "first-overlapping-group")
    gv.add_group_label(2, 5, "second-overlapping-group")
    with pytest.warns(RuntimeWarning, match="too narrow|fig_width|smaller labels"):
        fig, axs, renderer = _plot(gv, fig_width=1.4)
    left = fig.subplotpars.left
    fig.canvas.draw()
    assert fig.subplotpars.left == pytest.approx(left, abs=1e-6)
    _assert_auto_geometry(
        fig,
        axs,
        fig.canvas.get_renderer(),
        _group_artists(fig),
        require_clearance=False,
        check_y=False,
    )


def test_hidden_ytick_labels_do_not_consume_gutter_space():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(
        gv, ["", ""], values=[_values(100.0), _values(200.0)], show_range=True
    )
    gv.add_group_label(0, 1, "grp")
    fig, axs, renderer = _plot(gv)
    line_before = _group_artists(fig)[0]["line"].get_window_extent(renderer).x0
    for ax in axs:
        ax.tick_params(axis="y", labelleft=False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    line_after = _group_artists(fig)[0]["line"].get_window_extent(renderer).x0
    assert line_after > line_before + 1.0


def test_ungrouped_long_name_defines_shared_gutter():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(
        gv,
        ["very_long_ungrouped_track_name", "B", "C", "D"],
        heights=[0.4] * 4,
    )
    gv.add_group_label(2, 3, "grp")
    fig, axs, renderer = _plot(gv)
    name_bb = axs[0].yaxis.label.get_window_extent(renderer)
    line_bb = _group_artists(fig)[0]["line"].get_window_extent(renderer)
    assert line_bb.x1 <= name_bb.x0 + _EPS_PX
    _assert_auto_geometry(fig, axs, renderer, _group_artists(fig))


def test_zero_is_manual_placement():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B"])
    gv.add_group_label(0, 1, "zero", x=0.0, x_line_offset=0.015)
    fig, axs, renderer = _plot(gv)
    item = _group_artists(fig)[0]
    assert item["text"].get_position()[0] == pytest.approx(-0.015, abs=1e-6)
    assert item["line"].get_xdata()[0] == pytest.approx(0.015, abs=1e-6)


def test_automatic_placement_ignores_x_line_offset():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B"])
    gv.add_group_label(0, 1, "auto", x_line_offset=0.25)
    fig, axs, renderer = _plot(gv)
    item = _group_artists(fig)[0]
    assert item["config"].x is None
    assert item["line"].get_xdata()[0] != pytest.approx(0.25, abs=1e-3)
    _assert_auto_geometry(fig, axs, renderer, [item])


def test_mixed_auto_and_manual_keep_manual_formula():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C", "D"], heights=[0.4] * 4)
    gv.add_group_label(0, 1, "auto")
    gv.add_group_label(2, 3, "manual", x=0.02, x_line_offset=0.015)
    fig, axs, renderer = _plot(gv)
    items = _group_artists(fig)
    auto_item = next(item for item in items if item["config"].x is None)
    manual_item = next(item for item in items if item["config"].x is not None)
    _assert_auto_geometry(fig, axs, renderer, [auto_item])
    assert manual_item["text"].get_ha() == "center"
    assert manual_item["text"].get_position()[0] == pytest.approx(0.005, abs=1e-6)
    assert manual_item["line"].get_xdata()[0] == pytest.approx(0.035, abs=1e-6)


def test_axis_marks_top_margin_survives_group_label_redraw():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C"], heights=[0.5] * 3)
    gv.add_group_label(0, 2, "all")
    fig, axs, _ = _plot(gv)
    top_after_plot = fig.subplotpars.top
    left_after_plot = fig.subplotpars.left
    gv.set_axis_marks([10, 20, 30], labels=["rs1", "rs2", "rs3"], label_fontsize=11)
    fig.canvas.draw()
    top_after_marks = fig.subplotpars.top
    fig.canvas.draw()
    fig.canvas.draw()
    assert fig.subplotpars.top == pytest.approx(top_after_marks, abs=1e-6)
    assert fig.subplotpars.left == pytest.approx(left_after_plot, abs=1e-6)
    if top_after_marks < top_after_plot - 1e-6:
        assert fig.subplotpars.top < top_after_plot


def test_tight_layout_wrapper_matches_flag():
    from matplotlib.layout_engine import TightLayoutEngine

    from pygv.group_label_layout import GroupLabelLayoutEngine

    gv_tight = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv_tight, ["A", "B", "C"])
    gv_tight.add_group_label(0, 2, "all")
    fig_tight, axs_tight, renderer_tight = _plot(
        gv_tight, force_tight_layout=True, height_scale_factor=1.5
    )
    engine_tight = fig_tight.get_layout_engine()
    assert isinstance(engine_tight, GroupLabelLayoutEngine)
    assert isinstance(engine_tight._inner, TightLayoutEngine)
    _assert_auto_geometry(fig_tight, axs_tight, renderer_tight, _group_artists(fig_tight))

    gv_loose = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv_loose, ["A", "B", "C"])
    gv_loose.add_group_label(0, 2, "all")
    fig_loose, axs_loose, renderer_loose = _plot(gv_loose, force_tight_layout=False)
    engine_loose = fig_loose.get_layout_engine()
    assert isinstance(engine_loose, GroupLabelLayoutEngine)
    assert engine_loose._inner is None
    _assert_auto_geometry(fig_loose, axs_loose, renderer_loose, _group_artists(fig_loose))


def test_plot_twice_does_not_duplicate_artists_on_new_figure():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B", "C"])
    gv.add_group_label(0, 2, "all")
    fig1, _, _ = _plot(gv)
    n1 = len(_group_artists(fig1))
    fig2, axs2, renderer2 = _plot(gv)
    assert fig2 is not fig1
    assert len(_group_artists(fig2)) == n1
    _assert_auto_geometry(fig2, axs2, renderer2, _group_artists(fig2))


def test_empty_group_label_is_ignored():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B"])
    with pytest.warns(RuntimeWarning, match="empty group label"):
        gv.add_group_label(0, 1, "   ")
    gv.add_group_label(0, 1, "ok")
    fig, axs, renderer = _plot(gv)
    items = _group_artists(fig)
    assert len(items) == 1
    assert items[0]["config"].label == "ok"
    _assert_auto_geometry(fig, axs, renderer, items)


def test_out_of_range_group_index_is_ignored():
    gv = GenomeViewer(n_ticks=3)
    _add_dynamic_tracks(gv, ["A", "B"])
    with pytest.warns(RuntimeWarning, match="out of range"):
        gv.add_group_label(0, 5, "too-far")
    gv.add_group_label(0, 1, "ok")
    fig, axs, renderer = _plot(gv)
    items = _group_artists(fig)
    assert len(items) == 1
    _assert_auto_geometry(fig, axs, renderer, items)


def test_removed_track_does_not_crash_stale_group_label():
    gv = GenomeViewer(n_ticks=3)
    tracks = _add_dynamic_tracks(gv, ["A", "B", "C"])
    gv.add_group_label(0, 2, "all")
    gv.remove_track(tracks[-1])
    fig, axs, renderer = _plot(gv)
    assert fig.get_layout_engine() is None or not _group_artists(fig)

