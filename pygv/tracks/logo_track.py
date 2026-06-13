import numpy as np
import pandas as pd
import pyBigWig
from pyfaidx import Fasta
from typing import ClassVar
from pydantic import Field, field_validator
from pygv.tracks.logomaker.Logo import Logo
from pygv.utils import check_accessibility

from .track import NumericalTrack, NumericalTrackConfig


class LogoTrackConfig(NumericalTrackConfig):
    values: object = Field(
        default=None, description="Logo matrix as a DataFrame or 2D array."
    )
    color_scheme: object = Field(
        default=None, description="Color scheme used to render logo characters."
    )
    font_name: str = Field(
        default="sans", description="Font name used to render logo characters."
    )
    stack_order: str = Field(
        default="big_on_top",
        description="Order used to stack logo glyphs: 'big_on_top', 'small_on_top', or 'fixed'.",
    )
    center_values: bool = Field(
        default=False, description="Whether to center each logo matrix row around zero."
    )
    flip_below: bool = Field(
        default=True, description="Whether to flip characters drawn below the x-axis."
    )
    shade_below: float = Field(
        default=0.0,
        ge=0,
        le=1,
        description="Amount of shading applied to characters below the x-axis.",
    )
    fade_below: float = Field(
        default=0.0,
        ge=0,
        le=1,
        description="Amount of fading applied to characters below the x-axis.",
    )
    fade_probabilities: bool = Field(
        default=False,
        description="Whether character opacity should follow probability values.",
    )

    @field_validator("stack_order")
    def _validate_stack_order(cls, value):
        if value in {"big_on_top", "small_on_top", "fixed"}:
            return value
        raise ValueError("stack_order must be 'big_on_top', 'small_on_top', or 'fixed'")

    @field_validator("values")
    def _validate_values(cls, value):
        if value is None or isinstance(value, pd.DataFrame):
            return value
        if isinstance(value, np.ndarray) and len(value.shape) == 2 and value.shape[1] in {4, 20}:
            return value
        raise ValueError("values must be None, a DataFrame, or a 2D array with 4 or 20 columns")


class DynseqTrackConfig(LogoTrackConfig):
    is_nucleotide: bool = Field(
        default=True, description="Whether dynamic sequence values should be treated as nucleotides."
    )


class LogoTrack(NumericalTrack, LogoTrackConfig):
    """
    Sequence logo track. The implementation is built on the top of logomaker.
    After creating the track object, logo matrix should be assigned using :attr:`values`.

    Parameters
    ----------
    track : str
        Placeholder
    kwargs :
        color_scheme : (str, dict, or array with length 3)
            Specification of logo colors. Default is 'gray'. Can take a variety of
            forms.
             - (str) A built-in Logomaker color scheme in which the color of each
             character is determined that character's identity. Options are,
                 + For DNA/RNA: 'classic', 'grays', or 'base_paring'.
                 + For protein: 'hydrophobicity', 'chemistry', or 'charge'.
             - (str) A built-in matplotlib color name such as 'k' or 'tomato'
             - (list) An RGB array, i.e., 3 floats with values in the interval [0,1]
             - (dict) A dictionary that maps characters to colors, E.g.,
                {'A': 'blue',
                 'C': 'yellow',
                 'G': 'green',
                 'T': 'red'}
        font_name: (str)
            The character font to use when rendering the logo. For a list of
            valid font names, run logomaker.list_font_names().

        stack_order: (str)
            Must be 'big_on_top', 'small_on_top', or 'fixed'. If 'big_on_top',
            stack characters away from x-axis in order of increasing absolute value.
            If 'small_on_top', stack glyphs away from x-axis in order of
            decreasing absolute value. If 'fixed', stack glyphs from top to bottom
            in the order that characters appear in the data frame.

        center_values: (bool)
            If True, the stack of characters at each position will be centered
            around zero. This is accomplished by subtracting the mean value
            in each row of the matrix from each element in that row.

        flip_below: (bool)
            If True, characters below the x-axis (which correspond to negative
            values in the matrix) will be flipped upside down.

        shade_below: (float in [0,1])
            The amount of shading to use for characters drawn below the x-axis.
            Larger numbers correspond to more shading (i.e., darker characters).

        fade_below: (float in [0,1])
            The amount of fading to use for characters drawn below the x-axis.
            Larger numbers correspond to more fading (i.e., more transparent
            characters).

        fade_probabilities: (bool)
            If True, the characters in each stack will be assigned an alpha value
            equal to their height. This option only makes sense if df is a
            probability matrix. For additional customization, use
            Logo.fade_glyphs_in_probability_logo().

    Raises
    ------
    ValueError will be raised if the len the values property is not equal to the span of plotting region as defined as `end` - `start`

    Examples
    --------

    .. plot:: ../examples/plot_logo.py
    """

    _FIELD_PRIVATE_ATTRS: ClassVar[dict] = {
        **NumericalTrack._FIELD_PRIVATE_ATTRS,
        "values": "_values",
    }

    def _sync_private_field(self, name):
        if name == "values":
            object.__setattr__(
                self, "_values", self._normalize_values(self.__dict__[name])
            )
            return
        super()._sync_private_field(name)

    @staticmethod
    def _normalize_values(value):
        if value is None:
            return None
        if isinstance(value, np.ndarray):
            if value.shape[1] == 4:
                return pd.DataFrame(value, columns=["A", "C", "G", "T"])
            if value.shape[1] == 20:
                return pd.DataFrame(
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
            raise ValueError(
                "When providing an array as pwm, "
                "the columns of the array must be standard nucleotides (4) "
                "or amino acids (20) sorted alphabetically. Otherwise, please "
                "provide the matrix as a DataFrame."
            )
        if isinstance(value, pd.DataFrame):
            return value
        raise ValueError("values must be None, a numpy array, or a pandas DataFrame")

    def __init__(self, track: str = "", **kwargs):
        config = LogoTrackConfig(**kwargs)
        super(LogoTrack, self).__init__(**kwargs)

        self._values = None
        self.values = config.values
        self.color_scheme = config.color_scheme
        self.font_name = config.font_name
        self.stack_order = config.stack_order
        self.center_values = config.center_values
        self.flip_below = config.flip_below
        self.shade_below = config.shade_below
        self.fade_below = config.fade_below
        self.fade_probabilities = config.fade_probabilities

    @property
    def values(self):
        """A matrix specifying character heights and positions.
        Rows correspond to positions while columns correspond to characters.
        If providing a numpy array, columns must be standard nucleotides (4)/amino acids (20)
        sorted alphabetically. If providing a pd.DataFrame, column names must be single
        characters and row indices must be integers.

        Returns
        -------
        pd.DataFrame :
            shape: sequence_len, acgt (4)
        """
        return self._values

    @values.setter
    def values(self, value):
        if isinstance(value, np.ndarray):
            if value.shape[1] == 4:  # DNA/RNA
                self._values = pd.DataFrame(value, columns=["A", "C", "G", "T"])
            elif value.shape[1] == 20:  # peptide/protein
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

    def _get(self, chromosome, start, end):
        xvalues = np.arange(start, end, step=1)
        if len(xvalues) != len(self.values):
            raise ValueError(
                "The length of the region (end-start) is different from values' length."
            )
        return xvalues, self.values

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw track

        Parameters
        ----------
        chromosome : str
            placeholder
        start : int
            placeholder
        end : int
            placeholder
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(LogoTrack, self)._draw_track(
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


class DynseqTrack(NumericalTrack, DynseqTrackConfig):
    """
    Dynseq-flavor sequence logo track. The implementation is built on the top of logomaker.
    After creating the track object, logo matrix should be assigned using :attr:`values`.

    Parameters
    ----------
    track : str or list[str]
        track file(s)
    seq_fasta : str
        Genome fasta file
    is_nucleotide : bool
        Whether the track is a nucleotide track
    kwargs :
        color_scheme : (str, dict, or array with length 3)
            Specification of logo colors. Default is 'gray'. Can take a variety of
            forms.
             - (str) A built-in Logomaker color scheme in which the color of each
             character is determined that character's identity. Options are,
                 + For DNA/RNA: 'classic', 'grays', or 'base_paring'.
                 + For protein: 'hydrophobicity', 'chemistry', or 'charge'.
             - (str) A built-in matplotlib color name such as 'k' or 'tomato'
             - (list) An RGB array, i.e., 3 floats with values in the interval [0,1]
             - (dict) A dictionary that maps characters to colors, E.g.,
                {'A': 'blue',
                 'C': 'yellow',
                 'G': 'green',
                 'T': 'red'}
        font_name: (str)
            The character font to use when rendering the logo. For a list of
            valid font names, run logomaker.list_font_names().

        stack_order: (str)
            Must be 'big_on_top', 'small_on_top', or 'fixed'. If 'big_on_top',
            stack characters away from x-axis in order of increasing absolute value.
            If 'small_on_top', stack glyphs away from x-axis in order of
            decreasing absolute value. If 'fixed', stack glyphs from top to bottom
            in the order that characters appear in the data frame.

        center_values: (bool)
            If True, the stack of characters at each position will be centered
            around zero. This is accomplished by subtracting the mean value
            in each row of the matrix from each element in that row.

        flip_below: (bool)
            If True, characters below the x-axis (which correspond to negative
            values in the matrix) will be flipped upside down.

        shade_below: (float in [0,1])
            The amount of shading to use for characters drawn below the x-axis.
            Larger numbers correspond to more shading (i.e., darker characters).

        fade_below: (float in [0,1])
            The amount of fading to use for characters drawn below the x-axis.
            Larger numbers correspond to more fading (i.e., more transparent
            characters).

        fade_probabilities: (bool)
            If True, the characters in each stack will be assigned an alpha value
            equal to their height. This option only makes sense if df is a
            probability matrix. For additional customization, use
            Logo.fade_glyphs_in_probability_logo().

    Raises
    ------
    ValueError will be raised if the len the values property is not equal to the span of plotting region as defined as `end` - `start`

    Examples
    --------

    .. plot:: ../examples/plot_logo.py
    """

    def __init__(
        self, track: str = "", seq_fasta: str = "", is_nucleotide: bool = True, **kwargs
    ):
        config = DynseqTrackConfig(is_nucleotide=is_nucleotide, **kwargs)
        super(DynseqTrack, self).__init__(**kwargs)
        self.is_nucleotide = config.is_nucleotide

        self.bw = []
        if isinstance(track, str):
            check_accessibility(track, allow_remote=True)
            self.bw.append(pyBigWig.open(track))
        else:
            for sub_track in track:
                check_accessibility(sub_track, allow_remote=True)
                self.bw.append(pyBigWig.open(sub_track))
        check_accessibility(seq_fasta, allow_remote=False)

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
        self._values = None
        self.seq_fasta = Fasta(seq_fasta)

        self.color_scheme = config.color_scheme
        self.font_name = config.font_name
        self.stack_order = config.stack_order
        self.center_values = config.center_values
        self.flip_below = config.flip_below
        self.shade_below = config.shade_below
        self.fade_below = config.fade_below
        self.fade_probabilities = config.fade_probabilities

    def _get(self, chromosome, start, end):
        xvalues = np.arange(start, end, step=1)
        values = np.stack(
            [_bw.values(chromosome, start, end, numpy=True) for _bw in self.bw]
        ).mean(axis=0)
        values = self.data_transform(values)
        seq = self.seq_fasta[chromosome][start:end]
        mat = np.zeros((end - start, len(self._voc)))
        for i, s in enumerate(seq):
            try:
                mat[i, self._voc.index(s)] = 1.0
            except:
                pass
        self._values = pd.DataFrame(mat * values[:, None], columns=self._voc)

        return xvalues, self._values

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        """
        Draw track

        Parameters
        ----------
        chromosome : str
            placeholder
        start : int
            placeholder
        end : int
            placeholder
        ax : :class:`matplotlib.pyplot.Axes`
            matplotlib.pyplot.Axes for this track
        index : int
            The first subplot (track), index==0, will have its top border and xticks shown up
        kwargs :

        Returns
        -------

        """
        super(DynseqTrack, self)._draw_track(
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
