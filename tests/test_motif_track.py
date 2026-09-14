"""Tests for MotifTrack and motif search helpers."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from pydantic import ValidationError

from pygv.errors.DataIntegrity import ChromosomeNotInReference
from pygv.tracks.motif_track import (
    MINUS_COLOR,
    PLUS_COLOR,
    MotifTrack,
    find_overlapping_matches,
    is_iupac_pattern,
    is_nucleotide_regex,
    motif_to_regex,
    reverse_complement,
    search_motif_hits,
)
from pygv.viewer import GenomeViewer


def _write_fasta(path: Path, sequence: str, name: str = "chr1") -> str:
    fasta = path / "ref.fa"
    fasta.write_text(f">{name}\n{sequence}\n")
    return str(fasta)


def test_iupac_and_regex_classification():
    assert is_iupac_pattern("AAAAAARNR")
    assert is_iupac_pattern("acgtn")
    assert not is_iupac_pattern("TATA[AT]A")
    assert is_nucleotide_regex("TATA[AT]A[AT]A")
    assert is_nucleotide_regex("TATAAA+")
    assert not is_nucleotide_regex("AAAAAARNR")


def test_motif_to_regex_iupac_and_literal():
    assert motif_to_regex("ACCGCT") == "ACCGCT"
    assert motif_to_regex("AAAAAARNR") == "AAAAAA[AG][ACGTN][AG]"
    assert motif_to_regex("TATA[AT]A[AT]A") == "TATA[AT]A[AT]A"
    assert motif_to_regex("TATAAA+") == "TATAAA+"


def test_invalid_motif_raises():
    with pytest.raises(ValueError, match="invalid"):
        motif_to_regex("RATA[AT]")
    with pytest.raises(ValidationError):
        MotifTrack.model_validate({"track": "x.fa", "motif": "RATA[AT]"})
    with pytest.raises(ValidationError):
        MotifTrack.model_validate({"track": "x.fa", "motif": ""})


def test_plus_strand_exact_and_overlapping_hits():
    hits = search_motif_hits("AAA", "GAAAAAG", "chr1", 10, "+")
    assert [(h.start, h.end, h.match) for h in hits] == [
        (11, 14, "AAA"),
        (12, 15, "AAA"),
        (13, 16, "AAA"),
    ]


def test_minus_strand_maps_to_positive_coordinates():
    # Plus sequence CAT is ATG on the minus strand.
    hits = search_motif_hits("ATG", "GGCATCC", "chr1", 100, "-")
    assert len(hits) == 1
    hit = hits[0]
    assert (hit.start, hit.end) == (102, 105)
    assert hit.strand == "-"
    assert hit.match == "ATG"


def test_reverse_complement_iupac():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAAAARNR") == "YNYTTTTTT"


def test_regex_tata_box_variants(tmp_path: Path):
    seq = ("GC" * 20) + "TATAAAAA" + ("GC" * 4) + "TATATAAA" + ("GC" * 4) + "TATAAATA"
    fasta = _write_fasta(tmp_path, seq)
    track = MotifTrack(fasta, "TATA[AT]A[AT]A", strand="+")
    hits = track.find_hits("chr1", 0, len(seq))
    assert [h.match for h in hits] == ["TATAAAAA", "TATATAAA", "TATAAATA"]
    assert hits[0].start == 40


def test_iupac_pattern_hits(tmp_path: Path):
    seq = "GCGCAAAAAAGCGGCGC"
    fasta = _write_fasta(tmp_path, seq)
    track = MotifTrack(fasta, "AAAAAARNR", strand="+")
    hits = track.find_hits("chr1", 0, len(seq))
    assert len(hits) == 1
    assert hits[0].match == "AAAAAAGCG"
    assert (hits[0].start, hits[0].end) == (4, 13)


def test_pair_defaults_and_both_strands(tmp_path: Path):
    seq = "TATAAAAAGCGCTTTTTATA"
    fasta = _write_fasta(tmp_path, seq)
    plus, minus = MotifTrack.pair(fasta, "TATAAAAA", name="TATA")
    assert plus.strand == "+"
    assert minus.strand == "-"
    assert plus.color == PLUS_COLOR
    assert minus.color == MINUS_COLOR
    assert plus.name == "TATA"
    assert minus.name == "TATA Negative"
    assert [h.start for h in plus.find_hits("chr1", 0, len(seq))] == [0]
    minus_hits = minus.find_hits("chr1", 0, len(seq))
    assert len(minus_hits) == 1
    assert minus_hits[0].start == seq.find("TTTTTATA")

    both = MotifTrack(fasta, "TATAAAAA", strand="both", name="both")
    both_hits = both.find_hits("chr1", 0, len(seq))
    assert {h.strand for h in both_hits} == {"+", "-"}


def test_missing_chromosome_raises(tmp_path: Path):
    fasta = _write_fasta(tmp_path, "ACGTACGT")
    track = MotifTrack(fasta, "ACGT")
    with pytest.raises(ChromosomeNotInReference):
        track.find_hits("chr2", 0, 8)


def test_plot_motif_hits(tmp_path: Path):
    seq = ("GC" * 30) + "TATAAAAA" + ("GC" * 20) + "TTTTTATA" + ("GC" * 10)
    fasta = _write_fasta(tmp_path, seq)
    gv = GenomeViewer()
    gv.add_tracks(MotifTrack.pair(fasta, "TATAAAAA"))
    gv.add_track(MotifTrack(fasta, "AAAAAARNR", strand="+", show_name=True))
    axs = gv.plot("chr1", 0, len(seq))
    assert len(axs) == 3
    plt.close("all")


def test_unknown_kwarg_and_strand_validation(tmp_path: Path):
    fasta = _write_fasta(tmp_path, "ACGT")
    with pytest.raises(ValidationError):
        MotifTrack(fasta, "ACGT", strand="forward")
    with pytest.raises(ValidationError):
        MotifTrack(fasta, "ACGT", colr="blue")


def test_fasta_constructor_alias(tmp_path: Path):
    fasta = _write_fasta(tmp_path, "ACGTACGT")
    track = MotifTrack(fasta=fasta, motif="ACGT", strand="+")
    assert track.track == fasta
    assert [hit.start for hit in track.find_hits("chr1", 0, 8)] == [0, 4]


def test_zero_width_regex_matches_are_skipped():
    matches = list(find_overlapping_matches("A*", "GAG"))
    assert [(m.start(), m.end(), m.group()) for m in matches] == [(1, 2, "A")]
    hits = search_motif_hits("A*", "GAG", "chr1", 0, "+")
    assert [(h.start, h.end, h.match) for h in hits] == [(1, 2, "A")]


def test_regex_hits_that_cross_window_boundaries(tmp_path: Path):
    fasta = _write_fasta(tmp_path, "GCAAAAGC")
    track = MotifTrack(fasta, "A+", strand="+")
    hits = track.find_hits("chr1", 5, 6)
    assert [(h.start, h.end, h.match) for h in hits] == [
        (2, 6, "AAAA"),
        (3, 6, "AAA"),
        (4, 6, "AA"),
        (5, 6, "A"),
    ]


def test_both_strands_hits_are_sorted_by_coordinate(tmp_path: Path):
    fasta = _write_fasta(tmp_path, "TTTGGGGGGGGGGGGGGGGGGGGAAA")
    track = MotifTrack(fasta, "AAA", strand="both")
    hits = track.find_hits("chr1", 0, 26)
    assert [(h.start, h.end, h.strand) for h in hits] == [
        (0, 3, "-"),
        (23, 26, "+"),
    ]
    track._pre_plot_hook("chr1", 0, 26)
    occupied = [lane for lane in track._lane_registries if lane.features]
    assert len(occupied) == 1
    assert [h.strand for h in occupied[0].features] == ["-", "+"]


def test_adjacent_half_open_hits_share_a_lane(tmp_path: Path):
    fasta = _write_fasta(tmp_path, "AAAGGG")
    track = MotifTrack(fasta, "AAA", strand="+")
    hits = track.find_hits("chr1", 0, 6)
    assert [(h.start, h.end) for h in hits] == [(0, 3)]
    track = MotifTrack(fasta, "AAA|GGG", strand="+")
    hits = track.find_hits("chr1", 0, 6)
    assert [(h.start, h.end) for h in hits] == [(0, 3), (3, 6)]
    track._pre_plot_hook("chr1", 0, 6)
    occupied = [lane for lane in track._lane_registries if lane.features]
    assert len(occupied) == 1
    assert [(h.start, h.end) for h in occupied[0].features] == [(0, 3), (3, 6)]
