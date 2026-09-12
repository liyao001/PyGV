"""Motif-hit track modeled on IGV's Motif Finder.

See https://igv.org/doc/desktop/#UserGuide/tools/motif_finder/
"""

from __future__ import annotations

import re
from collections import namedtuple
from typing import Any, Iterable, Optional

import numpy as np
from matplotlib.patches import Rectangle
from pydantic import AliasChoices, Field, PrivateAttr, field_validator, model_validator
from pyfaidx import Fasta

from pygv.errors.DataIntegrity import ChromosomeNotInReference
from pygv.utils import check_accessibility

from .bed_track import _LaneRegistry
from .track import AnnotationTrack
from .types import Color, MotifStrand, ShowMode

# IUPAC nucleic-acid codes -> regex, matching IGV's motif finder.
# See http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
_IUPAC_TO_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[CG]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGTN]",
}

_IUPAC_CHARS = frozenset(_IUPAC_TO_REGEX)
_REGEX_BASES = frozenset("ACGTN")
_COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")

PLUS_COLOR = "#1F77B4"
MINUS_COLOR = "#D62728"
# Extra bases fetched on each side so variable-length regex hits that
# start before or extend past the window are still found.
_REGEX_BOUNDARY_PAD = 1024

MotifHit = namedtuple("MotifHit", ["contig", "start", "end", "name", "strand", "match"])


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a nucleotide string."""
    return sequence.translate(_COMPLEMENT)[::-1]


def is_iupac_pattern(pattern: str) -> bool:
    """Return True if every character is an IUPAC nucleotide code."""
    if not pattern:
        return False
    return all(char in _IUPAC_CHARS for char in pattern.upper())


def is_nucleotide_regex(pattern: str) -> bool:
    """Return True if ``pattern`` is a valid regex using only A/C/G/T/N letters."""
    try:
        re.compile(pattern)
    except re.error:
        return False
    return all((not char.isalpha()) or char.upper() in _REGEX_BASES for char in pattern)


def motif_to_regex(pattern: str) -> str:
    """Convert an IUPAC motif or nucleotide regex to a Python regex string.

    IUPAC patterns (e.g. ``AAAAAARNR``) are expanded to character classes.
    Nucleotide regular expressions (e.g. ``TATA[AT]A[AT]A`` or ``TATAAA+``)
    are returned unchanged aside from treating ``U`` as ``T``.
    """
    stripped = pattern.strip()
    if not stripped:
        raise ValueError("motif must be a non-empty nucleotide pattern")
    if is_iupac_pattern(stripped):
        return "".join(_IUPAC_TO_REGEX[char] for char in stripped.upper())
    if is_nucleotide_regex(stripped):
        return stripped.replace("U", "T").replace("u", "T")
    raise ValueError(
        "Invalid motif. Use a nucleotide sequence, IUPAC ambiguity codes "
        "(e.g. AAAAAARNR), or a nucleotide regex (e.g. TATA[AT]A[AT]A). "
        f"{pattern!r} is invalid."
    )


def _normalize_dna(sequence: str) -> str:
    return sequence.upper().replace("U", "T")


def find_overlapping_matches(pattern: str, sequence: str) -> Iterable[re.Match[str]]:
    """Yield overlapping regex matches, matching IGV's Matcher.find(last + 1).

    Zero-width matches (e.g. ``A*`` at a non-A base) are skipped so they
    cannot flood a track with empty intervals.
    """
    regex = re.compile(pattern, re.IGNORECASE)
    start = 0
    length = len(sequence)
    while start <= length:
        match = regex.search(sequence, start)
        if match is None:
            break
        if match.start() < match.end():
            yield match
        start = match.start() + 1


def search_motif_hits(
    pattern: str,
    sequence: str,
    chromosome: str,
    pos_start: int,
    strand: str,
) -> list[MotifHit]:
    """Search ``sequence`` (positive-strand bases starting at ``pos_start``).

    ``strand`` must be ``+`` or ``-``. Negative-strand search reverse-complements
    the sequence and maps coordinates back to the positive strand, as in IGV.
    """
    if strand not in {"+", "-"}:
        raise ValueError("strand must be '+' or '-'")
    dna = _normalize_dna(sequence)
    if strand == "-":
        dna = reverse_complement(dna)
    seq_len = len(dna)
    hits = []
    for match in find_overlapping_matches(pattern, dna):
        if strand == "+":
            start = pos_start + match.start()
            end = pos_start + match.end()
        else:
            start = pos_start + seq_len - match.end()
            end = pos_start + seq_len - match.start()
        matched = match.group(0).upper().replace("U", "T")
        hits.append(
            MotifHit(
                contig=chromosome,
                start=start,
                end=end,
                name=matched,
                strand=strand,
                match=matched,
            )
        )
    if strand == "-":
        hits.reverse()
    return hits


class MotifTrack(AnnotationTrack):
    """Display motif hits in a reference FASTA, like IGV's Motif Finder.

    The results are interval features. Search one strand, or both strands in a
    single track. :meth:`pair` builds the IGV-style plus (blue) and minus (red)
    track pair.

    Patterns may be:

    * a nucleotide sequence, e.g. ``ACCGCT``
    * IUPAC ambiguity codes, e.g. ``AAAAAARNR``
    * a regular expression of nucleotides, e.g. ``TATA[AT]A[AT]A`` or ``TATAAA+``

    Parameters
    ----------
    track
        Path to the reference FASTA (``fasta`` / ``seq_fasta`` are aliases).
    motif
        Nucleotide, IUPAC, or regex pattern to search.
    """

    track: str = Field(
        description="Path to the reference FASTA file",
        validation_alias=AliasChoices("track", "fasta", "seq_fasta"),
    )
    motif: str = Field(description="Nucleotide, IUPAC, or regex motif")
    strand: MotifStrand = Field(
        default="+",
        kw_only=True,
        description="Search the plus strand, minus strand, or both",
    )
    show_mode: ShowMode = Field(
        default="expanded",
        kw_only=True,
        description=(
            "Collapse overlapping hits (`collapsed`) or keep them on separate "
            "lanes (`expanded`), matching IGV's default expanded motif tracks."
        ),
    )
    color: Color = Field(
        default=None,
        kw_only=True,
        description="Fill color for plus-strand hits (and minus hits when strand='-')",
    )
    minus_color: Color = Field(
        default=MINUS_COLOR,
        kw_only=True,
        description="Fill color for minus-strand hits when strand='both'",
    )
    height: float = Field(
        default=0.8,
        gt=0,
        kw_only=True,
        description="Height of each feature lane",
    )
    show_name: bool = Field(
        default=False,
        kw_only=True,
        description="Label each hit with the matched sequence",
    )
    show_arrows: bool = Field(
        default=True,
        kw_only=True,
        description="Draw strand arrows on each hit",
    )

    _fasta: Any = PrivateAttr(default=None)
    _regex_pattern: str = PrivateAttr(default="")
    _iupac_length: Optional[int] = PrivateAttr(default=None)
    _small_relative: float = PrivateAttr(default=0)

    def __init__(
        self, track: str | None = None, motif: str | None = None, **data: Any
    ) -> None:
        if motif is not None:
            data["motif"] = motif
        if track is not None:
            data.setdefault("track", track)
        path = data.pop("track", None)
        if path is None:
            path = data.pop("fasta", None)
        if path is None:
            path = data.pop("seq_fasta", None)
        data.pop("fasta", None)
        data.pop("seq_fasta", None)
        super().__init__(track=path, **data)

    @field_validator("motif")
    @classmethod
    def _validate_motif(cls, value: str) -> str:
        motif_to_regex(value)
        return value.strip()

    @model_validator(mode="after")
    def _apply_defaults(self):
        self._regex_pattern = motif_to_regex(self.motif)
        self._iupac_length = (
            len(self.motif) if is_iupac_pattern(self.motif) else None
        )
        if self.color is None:
            self.color = MINUS_COLOR if self.strand == "-" else PLUS_COLOR
        if not self.name:
            self.name = (
                f"{self.motif} Negative" if self.strand == "-" else self.motif
            )
        return self

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        check_accessibility(self.track, allow_remote=False)
        self._fasta = Fasta(self.track)

    @classmethod
    def pair(
        cls,
        track: str,
        motif: str,
        *,
        name: Optional[str] = None,
        plus_name: Optional[str] = None,
        minus_name: Optional[str] = None,
        plus_color: Any = PLUS_COLOR,
        minus_color: Any = MINUS_COLOR,
        **kwargs: Any,
    ) -> tuple[MotifTrack, MotifTrack]:
        """Return plus- and minus-strand tracks, as IGV's Motif Finder does.

        Default colors are blue for the plus strand and red for the minus strand.
        """
        kwargs.pop("strand", None)
        kwargs.pop("color", None)
        label = name or motif
        plus = cls(
            track,
            motif,
            strand="+",
            name=plus_name or label,
            color=plus_color,
            **kwargs,
        )
        minus = cls(
            track,
            motif,
            strand="-",
            name=minus_name or f"{label} Negative",
            color=minus_color,
            **kwargs,
        )
        return plus, minus

    def _chrom_record(self, chromosome: str):
        if chromosome in self._fasta:
            return self._fasta[chromosome]
        raise ChromosomeNotInReference(
            f"{chromosome!r} is not in {self.track}. "
            f"Available: {', '.join(list(self._fasta.keys())[:12])}"
        )

    def find_hits(self, chromosome: str, start: int, end: int) -> list[MotifHit]:
        """Return motif hits overlapping ``[start, end)`` on the requested strand(s)."""
        record = self._chrom_record(chromosome)
        chrom_len = len(record)
        pad = (
            self._iupac_length - 1
            if self._iupac_length
            else _REGEX_BOUNDARY_PAD
        )
        fetch_start = max(0, int(start) - pad)
        fetch_end = min(chrom_len, int(end) + pad)
        if fetch_start >= fetch_end:
            return []
        sequence = str(record[fetch_start:fetch_end])
        strands = ("+", "-") if self.strand == "both" else (self.strand,)
        hits = []
        for strand in strands:
            hits.extend(
                search_motif_hits(
                    self._regex_pattern,
                    sequence,
                    chromosome,
                    fetch_start,
                    strand,
                )
            )
        hits = [hit for hit in hits if hit.end > start and hit.start < end]
        hits.sort(key=lambda hit: (hit.start, hit.end, hit.strand))
        return hits

    def _get(self, chromosome, start, end):
        yield from self.find_hits(chromosome, start, end)

    def _hit_color(self, strand: str):
        if self.strand == "both" and strand == "-":
            return self.minus_color
        return self.color

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        super()._pre_plot_hook(chromosome, start, end, **kwargs)
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
                key = (visible_start, visible_end, interval.strand)
                if key in added:
                    continue
                added.add(key)

            if len(self._lane_registries) == 0:
                self._lane_registries.append(_LaneRegistry())

            for lane in self._lane_registries:
                if lane.max_coord is None:
                    active_lane = lane.offset
                    lane.min_coord = start_loc
                    lane.max_coord = end_loc
                    lane.features.append(interval)
                    break
                if (
                    lane.max_coord <= start_loc - text_padding
                    or self.show_mode == "collapsed"
                ):
                    active_lane = lane.offset
                    lane.min_coord = min(lane.min_coord, start_loc)
                    lane.max_coord = max(lane.max_coord, end_loc)
                    lane.features.append(interval)
                    break
                active_lane = None

            if (
                isinstance(self.allowed_feature_lanes, int)
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

        if not self._lane_registries:
            self._lane_registries.append(_LaneRegistry())
        elif len(self._lane_registries) == 1:
            self._lane_registries.append(
                _LaneRegistry(
                    offset=len(self._lane_registries),
                    min_coord=start,
                    max_coord=end,
                    features=[],
                )
            )

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        self._ax.set_xlim((start, end))
        self._small_relative = 0.004 * (end - start)

        for lane in self._lane_registries:
            empty_lane = True
            real_active_line = (self.patch_height + self.lane_space) * lane.offset
            for interval in lane.features:
                empty_lane = False
                start_loc = int(interval.start)
                end_loc = int(interval.end)
                visible_start = max(start_loc, start)
                visible_end = min(end_loc, end)
                color = self._hit_color(interval.strand)

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
                rec = Rectangle(
                    xy=(start_loc, -1 * real_active_line - (self.patch_height / 2)),
                    width=max(end_loc - start_loc, 0),
                    height=self.patch_height,
                    edgecolor=self.edge_color,
                    facecolor=color,
                    linewidth=self.line_width,
                    clip_on=True,
                )
                self._ax.add_patch(rec)

                if self.show_arrows and interval.strand in {"+", "-"}:
                    step = int(self.arrow_interval * self._small_relative)
                    span = visible_end - visible_start
                    if span > 0:
                        if end_loc - start_loc > self._small_relative and step > 0:
                            positions = np.arange(
                                visible_start + self._small_relative,
                                visible_end + self._small_relative,
                                step,
                            )
                        else:
                            positions = (
                                (visible_start + visible_end) / 2,
                            )
                        for xpos in positions:
                            self._plot_gene_direction(
                                ax, xpos, -1 * real_active_line, interval.strand
                            )

                if self.show_name and interval.name:
                    if start_loc > start and interval.strand == "+":
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
                    elif end_loc < end and interval.strand == "-":
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
        self._ax.set_yticks([], minor=True)
        if index != 0:
            self._ax.set_xticks([])
            self._ax.set_xticks([], minor=True)
