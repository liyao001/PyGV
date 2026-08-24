"""Regression tests for v2 loader, logo, highlight, and layout fixes."""

from pathlib import Path

import numpy as np
import pandas as pd

from pygv.tracks.bed_track import BedPETrack
from pygv.tracks.gtf_track import GtfTrack
from pygv.tracks.logo_track import LogoTrack
from pygv.tracks.track import Track


def test_logo_ndarray_values_stay_dataframe():
    logo = LogoTrack("")
    matrix = np.array(
        [
            [0.25, 0.25, 0.25, 0.25],
            [0.1, 0.2, 0.3, 0.4],
        ]
    )
    logo.values = matrix
    assert isinstance(logo.values, pd.DataFrame)
    assert list(logo.values.columns) == ["A", "C", "G", "T"]
    assert np.allclose(logo.values.to_numpy(), matrix)


def test_add_highlight_region_keeps_color_alpha():
    track = Track(name="demo")
    track.add_highlight_region(10, 20)
    assert track._highlight_starts == [10.0]
    assert track._highlight_ends == [20.0]
    assert len(track._highlight_colors) == 1
    assert len(track._highlight_alphas) == 1
    track.remove_highlight()
    assert track._highlight_starts == []
    assert track._highlight_colors == []
    assert track._highlight_alphas == []


def test_uncompressed_gtf_uses_pandas_parser(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\ttranscript\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\ttest\texon\t160\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    track = GtfTrack(str(gtf), name="genes", height=0.5)
    assert track._parser == track._pd_parser
    records = list(track._get("chr1", 0, 300))
    assert len(records) == 1
    rec = records[0]
    assert rec.transcript_id == "t1"
    assert rec.start == 100
    assert rec.end == 200
    assert rec.exons == [(100, 150), (159, 200)]

    track._pre_plot_hook("chr1", 0, 300)
    assert track.height == 0.5
    assert track.layout_height() == 0.5


def test_uncompressed_bedpe_uses_pandas_parser(tmp_path: Path):
    bedpe = tmp_path / "loops.bedpe"
    bedpe.write_text("chr1\t10\t20\tchr1\t50\t60\tloop1\t.\t+\t-\n")
    track = BedPETrack(str(bedpe), name="loops")
    assert track._parser == track._pd_parser
    records = list(track._get("chr1", 0, 100))
    assert len(records) == 1
    assert records[0].name == "loop1"
