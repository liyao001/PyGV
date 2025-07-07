from typing import Union
from warnings import warn

import matplotlib as mpl
import matplotlib.cm
import matplotlib.font_manager
import numpy as np

import pygv.tracks
from pygv import __version__


class GenomeViewer(object):
    """
    Genome Viewer

    Examples
        >>> from pygv.viewer import GenomeViewer
        >>> from pygv.tracks import gtf_track
        >>> gv = GenomeViewer()
        >>> gencode_track = gtf_track.GtfTrack(
        >>>     "~/gencode.v34lift37.annotation.sorted.gtf.gz",
        >>>     name="GENCODE", show_genes=True, show_transcript_id=True,
        >>>     filters=lambda x: x.transcript_id in {"ENST00000332995.11_1", "ENSG00000112137.17_4",
        >>>                                           "ENST00000379350.5_1", "ENST00000379335.7_1"},
        >>>     annotation_formatter=lambda x: x.split(".")[0]
        >>> )
        >>> gv.add_track(gencode_track)
        >>> gv.plot("chr6", 12714999, 13292716)
        >>> plt.show()
    """

    def __init__(
        self,
        font_name=None,
        font_size=None,
        alternative_color_map=None,
        hspace=0.2,
        inward_ticks=None,
        n_ticks=None,
    ):
        """
        Initiate a new genome viewer instance

        Parameters
        ----------
        font_name : None or str
             Name of a font, which should be included in `matplotlib.font_manager.fontManager.ttflist`
        font_size : None or int
             Default size for texts
        alternative_color_map : None or str

        hspace : float
            The amount of height reserved for space between tracks,
            expressed as a fraction of the average axis height. By default, 0.2.
        inward_ticks : bool or None
            Set this as True or False if you want to override the `inward_ticks` properties of
            each individual track.
        n_ticks : int or None
            Number of ticks for the x-axis. If None, number of ticks will be determined automatically.

        """
        self._registered_tracks = []
        if font_name is not None and font_name in GenomeViewer._supported_fonts():
            mpl.pyplot.rcParams["font.family"] = font_name

        if font_size is not None:
            mpl.pyplot.rcParams["font.size"] = font_size

        self.alternative_colors = []
        if alternative_color_map is not None and type(alternative_color_map) is str:
            self.alternative_colors = matplotlib.cm.get_cmap(
                alternative_color_map
            ).colors

        self._plot_chrom = None
        self._plot_start = None
        self._plot_end = None
        self._hspace = hspace
        self._inward_ticks = inward_ticks
        self._n_ticks = n_ticks
        self._group_auto_scales = []

    @staticmethod
    def _supported_fonts():
        try:
            return sorted(
                set([f.name for f in matplotlib.font_manager.fontManager.ttflist])
            )
        except AttributeError:
            return sorted(
                set([f._name for f in matplotlib.font_manager.fontManager.ttflist])
            )

    def add_track(self, track: pygv.tracks.track.Track) -> None:
        """
        Add a track to a `GenomeViewer` instance

        Parameters
        ----------
        track : tracks.track.Track
            Track object to be added
        Returns
        -------

        """
        self._registered_tracks.append(track)

    def add_group_autoscale(self, track_idx: Union[tuple[int, ...], list[int]]):
        """
        Add group autoscale

        Parameters
        ----------
        track_idx : Union[tuple[int, ...], list[int]]
            Indexes of the tracks to be scaled together

        Examples
        --------

        .. plot:: ../examples/plot_group_autoscale.py
        """
        tracks = []
        for tid in track_idx:
            if tid < len(self._registered_tracks):
                if isinstance(self._registered_tracks[tid], pygv.tracks.NumericalTrack):
                    tracks.append(tid)
                else:
                    warn(f"Track {tid} is not scalable", RuntimeWarning)
            else:
                warn(
                    f"Track index {tid} is larger than the number of registered tracks",
                    RuntimeWarning,
                )
        if len(tracks) > 0:
            self._group_auto_scales.append(tracks)

    def add_group_autoscale_by_name(
        self, track_name: Union[tuple[str, ...], list[str]]
    ):
        """
        Add group autoscale

        Parameters
        ----------
        track_name : Union[tuple[str, ...], list[str]]
            Names of the tracks to be scaled together

        """
        tracks = []
        all_track_names = [t.name for t in self._registered_tracks]
        for tname in track_name:
            try:
                tid = all_track_names.index(tname)
                if isinstance(self._registered_tracks[tid], pygv.tracks.NumericalTrack):
                    tracks.append(tid)
                else:
                    warn(f"Track {tid} is not scalable", RuntimeWarning)
            except ValueError:
                warn(f"Cannot find Track {tname}", RuntimeWarning)

        if len(tracks) > 0:
            self._group_auto_scales.append(tracks)

    def add_tracks(self, tracks):
        """
        Add tracks to a `GenomeViewer` instance

        Parameters
        ----------
        tracks : tuple or list
            Objects of `tracks.track.Track` to be added
        Returns
        -------

        """
        for track in tracks:
            self._registered_tracks.append(track)

    def remove_track(self, track):
        """
        Remove a track from a `GenomeViewer` instance

        Parameters
        ----------
        track : tracks.track.Track
            Track object to be removed

        Returns
        -------

        """
        if track in self._registered_tracks:
            self._registered_tracks.remove(track)

    def reset_group_autoscale(self):
        """
        Remove all group autoscale rules
        """
        self._group_auto_scales = []

    def set_highlight_regions(
        self,
        starts: Union[list, tuple],
        ends: Union[list, tuple],
        colors=(),
        alpha_vals=(),
    ):
        """
        Set highlight regions for all tracks. If you only want to highlight regions on specific tracks,
        you can call each track's :meth:`~pygv.tracks.track.Track.set_highlight_regions` method.
        Chromosome name is not needed for this method, it will use the same chromosome name when
        you call the :meth:`~pygv.viewer.GenomeViewer.plot` method.

        Parameters
        ----------
        starts : Union[list, tuple]
            Start positions
        ends : Union[list, tuple]
            End positions
        colors : tuple
            Leave it as an empty tuple if you want to use the default color.
            If you only give one color, it will be applied to all regions; otherwise, you should specify
            colors for each region.
        alpha_vals : tuple
            Leave it as an empty tuple if you want to use the default transparency level (0.5).
            If you only give one value, it will be applied to all regions; otherwise, you should specify
            transparency values for each region.

        Returns
        -------

        Examples
        --------

        .. plot:: ../examples/plot_highlight_regions.py
        """
        if len(self._registered_tracks) == 0:
            raise RuntimeWarning(
                "You need to add tracks before adding highlight regions"
            )

        for track in self._registered_tracks:
            track.set_highlight_regions(starts, ends, colors, alpha_vals)

    def set_global_highlight_region(
        self, start: int, end: int, color="yellow", alpha=0.3
    ):
        """
        Set a global highlight region across all tracks, including the spaces between subplots.

        Parameters
        ----------
        start : int
            Start position of the highlight region (genomic coordinate).
        end : int
            End position of the highlight region (genomic coordinate).
        color : str, optional
            Color of the highlight region. Default is "yellow".
        alpha : float, optional
            Transparency level of the highlight region. Default is 0.3.

        Returns
        -------
        None
        """
        if (
            self._plot_chrom is None
            or self._plot_start is None
            or self._plot_end is None
        ):
            raise RuntimeError(
                "You must call the `plot` method before adding a global highlight region."
            )

        # Normalize the start and end positions to the x-axis range
        x_start = max(self._plot_start, start)
        x_end = min(self._plot_end, end)

        if x_start >= x_end:
            warn(
                "The highlight region is outside the plotted range and will not be displayed.",
                RuntimeWarning,
            )
            return

        # Add a rectangle to the background of the figure
        fig = mpl.pyplot.gcf()
        ax = fig.add_subplot(111, label="global_highlight", zorder=-1, frame_on=False)
        ax.set_xlim(self._plot_start, self._plot_end)
        ax.set_ylim(0, 1)
        ax.axis("off")  # Hide axes

        # Add the rectangle
        rect = mpl.patches.Rectangle(
            (x_start, 0),
            x_end - x_start,
            1,
            color=color,
            alpha=alpha,
            transform=ax.transData,
            zorder=-1,
        )
        ax.add_patch(rect)

    def show_tracks(self):
        """
        Show all registered tracks

        Returns
        -------
        tracks : list
            A list of registered tracks. Each element is also a list: name of the track, track type, track.
        """
        tracks = []
        for track in self._registered_tracks:
            tracks.append((track.name, type(track), track))
        return tracks

    def plot(
        self,
        chromosome,
        start,
        end,
        fig_width=8,
        height_scale_factor=1,
        force_tight_layout=None,
        fig_height=None,
        **kwargs,
    ):
        """
        Plot the genome viewer with the registered tracks.

        Parameters
        ----------
        chromosome : str
            Chromosome/contig the region locates.
        start : int
            Start of the genomic region, 0-based.
        end : int
            End of the genomic region, 0-based.
        fig_width : float, optional
            Width (in inches) of the figure. Default is 8.
        height_scale_factor : float, optional
            Aspect ratio of the figure, so that ``height_scale_factor`` * ``fig_width`` gives the height of the figure.
            Default is 1.
        force_tight_layout : bool or None, optional
            If True, PyGV applies tight layout to the figure. Default is None.
        fig_height : float or None, optional
            Height of the figure. If None, height will be the sum of tracks' heights (in unit) * ``height_scale_factor``.
            Default is None.
        **kwargs : dict, optional
            Additional keyword arguments for track customization.

        Returns
        -------
        list of matplotlib.pyplot.Axes
            Axes for each track.
        """
        self._plot_chrom = chromosome
        self._plot_start = start
        self._plot_end = end

        # Validate and prepare tracks
        self._validate_tracks()
        self._prepare_tracks(chromosome, start, end)

        # Calculate figure dimensions
        heights = self._calculate_track_heights()
        fig_height = self._determine_figure_height(
            heights, height_scale_factor, fig_height
        )

        # Create figure and axes
        fig, axs = self._create_figure_and_axes(fig_width, fig_height, heights)

        # Draw tracks
        self._draw_tracks(axs, chromosome, start, end, **kwargs)

        # Adjust layout and apply group autoscale
        self._adjust_layout(fig, axs, force_tight_layout)

        return axs

    def _validate_tracks(self):
        """
        Ensure that tracks are registered before plotting.

        Raises
        ------
        RuntimeError
            If no tracks are registered.
        """
        if len(self._registered_tracks) == 0:
            raise RuntimeError(
                "No tracks registered, please add tracks to the Viewer first."
            )

    def _prepare_tracks(self, chromosome, start, end):
        """
        Prepare tracks for plotting by calling their pre-plot hooks.

        Parameters
        ----------
        chromosome : str
            Chromosome/contig the region locates.
        start : int
            Start of the genomic region, 0-based.
        end : int
            End of the genomic region, 0-based.
        """
        for track in self._registered_tracks:
            track._pre_plot_hook(
                chromosome=chromosome,
                start=start,
                end=end,
                inward_ticks=self._inward_ticks,
            )

    def _calculate_track_heights(self):
        """
        Calculate the heights of all registered tracks.

        Returns
        -------
        numpy.ndarray
            Array of track heights.
        """
        heights = np.zeros(len(self._registered_tracks))
        for index, track in enumerate(self._registered_tracks):
            heights[index] = track.height
        return heights

    def _determine_figure_height(self, heights, height_scale_factor, fig_height):
        """
        Determine the height of the figure.

        Parameters
        ----------
        heights : numpy.ndarray
            Array of track heights.
        height_scale_factor : float
            Aspect ratio of the figure.
        fig_height : float or None
            Predefined figure height.

        Returns
        -------
        float
            Calculated figure height.
        """
        if fig_height is None:
            fig_height = heights.sum() * height_scale_factor
        return fig_height

    def _create_figure_and_axes(self, fig_width, fig_height, heights):
        """
        Create the matplotlib figure and axes.

        Parameters
        ----------
        fig_width : float
            Width of the figure.
        fig_height : float
            Height of the figure.
        heights : numpy.ndarray
            Array of track heights.

        Returns
        -------
        tuple
            A tuple containing the figure and a list of axes.
        """
        normed_heights = heights
        fig, axs = mpl.pyplot.subplots(
            figsize=(fig_width, fig_height),
            ncols=1,
            nrows=len(self._registered_tracks),
            gridspec_kw={"height_ratios": normed_heights},
        )
        if isinstance(axs, mpl.pyplot.Axes):
            axs = [axs]
        return fig, axs

    def _draw_tracks(self, axs, chromosome, start, end, **kwargs):
        """
        Draw each track on the corresponding axis.

        Parameters
        ----------
        axs : list of matplotlib.pyplot.Axes
            Axes for each track.
        chromosome : str
            Chromosome/contig the region locates.
        start : int
            Start of the genomic region, 0-based.
        end : int
            End of the genomic region, 0-based.
        **kwargs : dict, optional
            Additional keyword arguments for track customization.
        """
        for index, track in enumerate(self._registered_tracks):
            sax = axs[index]
            track._draw_track(
                chromosome=chromosome,
                start=start,
                end=end,
                ax=sax,
                index=index,
                n_ticks=self._n_ticks,
                **kwargs,
            )
            track._post_plot_hook(chromosome, start, end, ax=sax, index=index, **kwargs)

    def _adjust_layout(self, fig, axs, force_tight_layout):
        """
        Adjust the layout and apply group autoscale.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            The figure object.
        axs : list of matplotlib.pyplot.Axes
            Axes for each track.
        force_tight_layout : bool or None
            Whether to apply tight layout.
        """
        mpl.pyplot.subplots_adjust(hspace=self._hspace)

        # Apply group autoscale
        if len(self._group_auto_scales) > 0:
            self._apply_group_autoscale()

        fig.align_ylabels()
        if force_tight_layout is None or not force_tight_layout:
            fig.set_tight_layout(False)
        else:
            fig.set_tight_layout(True)

    def _apply_group_autoscale(self):
        """
        Apply group autoscale to the tracks.
        """
        for group in self._group_auto_scales:
            y_lims = [self._registered_tracks[t]._ax.get_ylim() for t in group]
            spans = [lim[1] - lim[0] for lim in y_lims]
            max_idx = np.argmax(spans)
            target_ylim = y_lims[max_idx]
            track_id = group[max_idx]
            target_yticks = self._registered_tracks[track_id]._ax.get_yticks()
            target_yticklabels = self._registered_tracks[track_id]._ax.get_yticklabels()
            target_yscale = self._registered_tracks[track_id]._ax.get_yscale()
            target_yscale_func = self._registered_tracks[track_id]._yscale_func
            for _track in group:
                if target_yscale == "function":
                    self._registered_tracks[_track]._ax.set_yscale(
                        target_yscale, functions=target_yscale_func
                    )
                self._registered_tracks[_track]._ax.set_ylim(target_ylim)
                self._registered_tracks[_track]._ax.set_yticks(target_yticks)
                self._registered_tracks[_track]._ax.set_yticklabels(target_yticklabels)

    def save(self, *args, **kwargs):
        """
        Save figure to a file

        Parameters
        ----------
        args
        kwargs

        Returns
        -------

        """
        import datetime
        import getpass

        metadata = None
        if args[0].find(".pdf") != -1:
            metadata = {
                "Title": "Genome viewer shot at {0}:{1}-{2}".format(
                    self._plot_chrom, self._plot_start, self._plot_end
                ),
                "Author": getpass.getuser(),
                "Creator": "Python Genome Viewer ver{0}".format(__version__),
                "Producer": "PyGV ver{0} via matplotlib ver{1}".format(
                    __version__, mpl.__version__
                ),
                "CreationDate": datetime.datetime.now(),
            }
            kwargs["metadata"] = metadata
            kwargs["bbox_inches"] = "tight"
        mpl.pyplot.savefig(*args, **kwargs)
