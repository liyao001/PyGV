from __future__ import annotations

from typing import Any, List, Union

import numpy as np
import pandas as pd
import pyBigWig
from pydantic import Field, PrivateAttr
from pyfaidx import Fasta

from pygv.tracks.logomaker.Logo import Logo
from pygv.utils import check_accessibility

from .track import NumericalTrack
from .types import StackOrder


class LogoTrack(NumericalTrack):
    """Sequence logo track built on logomaker. Assign the logo matrix via `values`."""

    track: str = Field(default="", description="Placeholder")
    color_scheme: Any = Field(
        default=None,
        kw_only=True,
        description="Logo colors: a logomaker scheme, matplotlib color, RGB array, or dict",
    )
    font_name: str = Field(default="sans", kw_only=True)
    stack_order: StackOrder = Field(default="big_on_top", kw_only=True)
    center_values: bool = Field(default=False, kw_only=True)
    flip_below: bool = Field(default=True, kw_only=True)
    shade_below: float = Field(default=0.0, ge=0, le=1, kw_only=True)
    fade_below: float = Field(default=0.0, ge=0, le=1, kw_only=True)
    fade_probabilities: bool = Field(default=False, kw_only=True)

    _values: Any = PrivateAttr(default=None)

    def __init__(self, track: str = "", **data: Any) -> None:
        super().__init__(track=track, **data)

    @property
    def values(self):
        """A matrix specifying character heights and positions."""
        return self._values

    @values.setter
    def values(self, value):
        if isinstance(value, np.ndarray):
            if value.shape[1] == 4:
                self._values = pd.DataFrame(value, columns=["A", "C", "G", "T"])
            elif value.shape[1] == 20:
                self._values = pd.DataFrame(
                    value,
                    columns=[
                        "A",
                        "C",
                        "D",
                        "E",
                        "F",
                        "G",
                        "H",
                        "I",
                        "K",
                        "L",
                        "M",
                        "N",
                        "P",
                        "Q",
                        "R",
                        "S",
                        "T",
                        "V",
                        "W",
                        "Y",
                    ],
                )
            else:
                raise ValueError(
                    "When providing an array as pwm, "
                    "the columns of the array must be standard nucleotides (4) "
                    "or amino acids (20) sorted alphabetically. Otherwise, please "
                    "provide the matrix as a DataFrame."
                )
        elif isinstance(value, pd.DataFrame):
            self._values = value
        else:
            raise TypeError(
                "values must be a numpy ndarray or pandas DataFrame, "
                f"got {type(value).__name__}"
            )

    def _get(self, chromosome, start, end):
        xvalues = np.arange(start, end, step=1)
        if len(xvalues) != len(self.values):
            raise ValueError(
                "The length of the region (end-start) is different from values' length."
            )
        return xvalues, self.values

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        x, df = self._get(chromosome, start, end)
        Logo(
            df,
            offset=start,
            ax=ax,
            color_scheme=self.color_scheme,
            font_name=self.font_name,
            stack_order=self.stack_order,
            center_values=self.center_values,
            flip_below=self.flip_below,
            shade_below=self.shade_below,
            fade_below=self.fade_below,
            fade_probabilities=self.fade_probabilities,
            baseline_width=self.line_width,
        )
        self._ax = ax


class DynseqTrack(NumericalTrack):
    """Dynseq-flavor sequence logo track. Assign signal via bigWig and sequence via FASTA."""

    track: Union[str, List[str]] = Field(
        default="", description="BigWig file path(s) used as letter heights"
    )
    seq_fasta: str = Field(default="", description="Genome FASTA file")
    is_nucleotide: bool = Field(default=True, description="Whether the alphabet is DNA")
    color_scheme: Any = Field(default=None, kw_only=True)
    font_name: str = Field(default="sans", kw_only=True)
    stack_order: StackOrder = Field(default="big_on_top", kw_only=True)
    center_values: bool = Field(default=False, kw_only=True)
    flip_below: bool = Field(default=True, kw_only=True)
    shade_below: float = Field(default=0.0, ge=0, le=1, kw_only=True)
    fade_below: float = Field(default=0.0, ge=0, le=1, kw_only=True)
    fade_probabilities: bool = Field(default=False, kw_only=True)

    _bw: list = PrivateAttr(default_factory=list)
    _voc: tuple = PrivateAttr(default=())
    _values: Any = PrivateAttr(default=None)
    _seq_fasta: Any = PrivateAttr(default=None)

    def __init__(
        self,
        track: Union[str, List[str]] = "",
        seq_fasta: str = "",
        is_nucleotide: bool = True,
        **data: Any,
    ) -> None:
        super().__init__(
            track=track, seq_fasta=seq_fasta, is_nucleotide=is_nucleotide, **data
        )

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        opened = []
        if isinstance(self.track, str):
            if self.track:
                check_accessibility(self.track, allow_remote=True)
                opened.append(pyBigWig.open(self.track))
        else:
            for sub_track in self.track:
                check_accessibility(sub_track, allow_remote=True)
                opened.append(pyBigWig.open(sub_track))
        self._bw = opened
        if self.seq_fasta:
            check_accessibility(self.seq_fasta, allow_remote=False)
            self._seq_fasta = Fasta(self.seq_fasta)
        if self.is_nucleotide:
            self._voc = ("A", "C", "G", "T")
        else:
            self._voc = (
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y",
            )

    def _get(self, chromosome, start, end):
        xvalues = np.arange(start, end, step=1)
        values = np.stack(
            [_bw.values(chromosome, start, end, numpy=True) for _bw in self._bw]
        ).mean(axis=0)
        values = self.data_transform(values)
        seq = self._seq_fasta[chromosome][start:end]
        mat = np.zeros((end - start, len(self._voc)))
        for i, s in enumerate(seq):
            try:
                mat[i, self._voc.index(s)] = 1.0
            except Exception:
                pass
        self._values = pd.DataFrame(mat * values[:, None], columns=self._voc)

        return xvalues, self._values

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super()._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        x, df = self._get(chromosome, start, end)
        Logo(
            df,
            offset=start,
            ax=ax,
            color_scheme=self.color_scheme,
            font_name=self.font_name,
            stack_order=self.stack_order,
            center_values=self.center_values,
            flip_below=self.flip_below,
            shade_below=self.shade_below,
            fade_below=self.fade_below,
            fade_probabilities=self.fade_probabilities,
            baseline_width=self.line_width,
        )
        self._ax = ax
