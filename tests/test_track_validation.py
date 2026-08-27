"""Validation tests for pydantic Track parameters."""

import pytest
from pydantic import ValidationError

from pygv.tracks.bed_track import BedTrack
from pygv.tracks.logo_track import LogoTrack
from pygv.tracks.track import (
    AnnotationTrack,
    DualAxisTrack,
    DynamicValueTrack,
    NumericalTrack,
    Track,
)


def test_track_default_constructs():
    track = Track(name="demo")
    assert track.name == "demo"
    assert track.height == 1
    assert track.alpha == 0.8
    assert track.layout_height() == 1


def test_invalid_alpha_raises():
    with pytest.raises(ValidationError):
        Track(alpha=1.5)
    track = Track()
    with pytest.raises(ValidationError):
        track.alpha = -0.1


def test_invalid_color_raises():
    with pytest.raises(ValidationError):
        Track(color="not-a-color")


def test_invalid_height_raises():
    with pytest.raises(ValidationError):
        Track(height=0)
    with pytest.raises(ValidationError):
        Track(height=-1)


def test_unknown_kwarg_raises():
    with pytest.raises(ValidationError):
        Track(colr="#000000")


def test_inward_ticks_alias():
    track = Track(inward_ticks=True)
    assert track.inward_yticks is True
    track = Track(inward_yticks=False)
    assert track.inward_yticks is False


def test_annotation_track_positional_path():
    track = AnnotationTrack("genes.gtf", name="genes", show_name=False)
    assert track.track == "genes.gtf"
    assert track.show_name is False
    assert track.layout_height() == 1


def test_numerical_transformation_alias():
    import numpy as np

    track = NumericalTrack()
    assert np.allclose(track.data_transform([1.0, 2.0]), [1.0, 2.0])
    track = NumericalTrack(data_transform=None)
    assert np.allclose(track.data_transform([1.0, 2.0]), [1.0, 2.0])
    track = NumericalTrack(transformation="asinh")
    assert np.allclose(track.data_transform([0.0]), np.arcsinh([0.0]))
    track.data_transform = "log2"
    assert np.allclose(track.data_transform([2.0]), np.log2([2.0]))


def test_invalid_show_mode_raises():
    with pytest.raises(ValidationError):
        BedTrack.model_validate({"track": "x.bed", "show_mode": "hidden"})


def test_logo_and_dynamic_tracks_construct():
    logo = LogoTrack("", name="logo", stack_order="small_on_top")
    assert logo.stack_order == "small_on_top"
    dyn = DynamicValueTrack("", name="Demo")
    assert dyn.name == "Demo"


def test_invalid_stat_and_transform_raise():
    with pytest.raises(ValidationError):
        NumericalTrack(stat_method="nope")
    with pytest.raises(ValidationError):
        NumericalTrack(transformation="not-a-transform")
    with pytest.raises(ValidationError):
        NumericalTrack(data_transform="not-a-transform")


def test_dual_axis_layout_height():
    left = Track(name="left", height=2)
    right = Track(name="right", height=3)
    dual = DualAxisTrack(left, right)
    assert dual.layout_height() == 3
    assert dual.left_track is left
    assert dual.right_track is right


def test_paired_strandless_accepts_pos_neg_color(monkeypatch):
    from pygv.tracks import bigwig_track

    class _FakeBw:
        def values(self, *args, **kwargs):
            return None

    monkeypatch.setattr(bigwig_track, "check_accessibility", lambda *a, **k: None)
    monkeypatch.setattr(bigwig_track.pyBigWig, "open", lambda path: _FakeBw())

    track = bigwig_track.PairedStrandlessTrack(
        "pl.bw",
        "mn.bw",
        pos_color="#ee6352",
        neg_color="#0077b6",
    )
    assert track.pos_color == "#ee6352"
    assert track.neg_color == "#0077b6"
    assert track.color == "#ee6352"
