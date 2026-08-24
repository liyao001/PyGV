from typing import Union
from warnings import warn

import matplotlib as mpl
import matplotlib.cm
import matplotlib.font_manager
import numpy as np

import pygv.tracks
from pygv import __version__
from pygv.configs.label import GroupLabelConfig, GroupLabels


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
        self._group_labels = GroupLabels()

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

    def add_group_label(self, start_track_idx: int, end_track_idx: int, label: str, x=0.02, x_line_offset=0.015):
        """
        Add group label

        Parameters
        ----------
        start_track_idx : int
            Index of the start track (0-based)
        end_track_idx : int
            Index of the end track (0-based)
        label : str
            Group label
        x : float
            X-position for the label in figure coordinates (default 0.02).
        x_line_offset : float
            Offset for the line from the label in figure coordinates (default 0.015).

        Examples
        --------

        .. plot:: ../examples/plot_group_label.py
        """
        n_total_tracks = len(self._registered_tracks)
        if end_track_idx >= n_total_tracks:
            warn(f"End track index {end_track_idx} out of range", RuntimeWarning)
        else:
            self._group_labels.add(GroupLabelConfig(
                start_track_idx=start_track_idx,
                end_track_idx=end_track_idx,
                label=label,
                x=x,
                x_line_offset=x_line_offset
            ))

    def add_group_label_by_name(
        self, start_track_name: str, end_track_name: str, label: str, x=0.02, x_line_offset=0.015
    ):
        """
        Add group label by track names

        Parameters
        ----------
        start_track_name : str
            Name of the start track
        end_track_name : str
            Name of the end track
        label : str
            Group label
        x : float
            X-position for the label in figure coordinates (default 0.02).
        x_line_offset : float
            Offset for the line from the label in figure coordinates (default 0.015).

        """
        all_track_names = [t.name for t in self._registered_tracks]
        try:
            stid = all_track_names.index(start_track_name)
            etid = all_track_names.index(end_track_name)
            self.add_group_label(stid, etid, label, x=x, x_line_offset=x_line_offset)
        except ValueError:
            warn(f"Cannot find Track(s) {start_track_name}, {end_track_name}", RuntimeWarning)

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

    def set_global_vertical_line(
        self,
        position: int,
        color="red",
        alpha=0.8,
        line_width=1.5,
        line_style="-",
        margin_frac=0.02,
    ):
        """
        Draw a global vertical line at a single genomic position across all tracks,
        including spaces between subplots.

        Parameters
        ----------
        position : int
            Genomic coordinate where the line should be drawn.
        color : str, optional
            Line color. Default is "red".
        alpha : float, optional
            Line transparency. Default is 0.8.
        line_width : float, optional
            Line width. Default is 1.5.
        line_style : str, optional
            Matplotlib line style. Default is "-".
        margin_frac : float, optional
            Fractional amount to extend the line beyond the top of the panel.
            Helps cover tiny plotting margins near the coordinate axis.
            Default is 0.02.
        """
        if (
            self._plot_chrom is None
            or self._plot_start is None
            or self._plot_end is None
        ):
            raise RuntimeError(
                "You must call the `plot` method before adding a global vertical line."
            )

        if not (self._plot_start <= position <= self._plot_end):
            warn(
                "The position is outside the plotted range and the vertical line will not be displayed.",
                RuntimeWarning,
            )
            return

        fig = mpl.pyplot.gcf()
        ax = fig.add_subplot(
            111,
            label=f"global_vline_{position}_{len(fig.axes)}",
            zorder=200,
            frame_on=False,
        )
        ax.set_xlim(self._plot_start, self._plot_end)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.axvline(
            x=position,
            ymin=0,
            ymax=1 + margin_frac,
            color=color,
            alpha=alpha,
            linewidth=line_width,
            linestyle=line_style,
            clip_on=False,
        )

    def set_axis_marks(
        self,
        positions,
        color="red",
        size=8,
        line_width=1.0,
        stem_length=0.06,
        labels=(),
        label_rotation=90,
        label_fontsize=None,
        label_color=None,
        label_offset=0.04,
    ):
        """
        Add SNP-style lollipop marks to the top coordinate spine of the first track.
        Each mark is a downward-pointing triangle (``"v"``) sitting on the spine with a
        short vertical stem dropping from it, matching the style used by the
        ``UCSCMutationTrack``. Optionally, a text label can be rendered centered above
        each mark.
        Must be called after :meth:`~pygv.viewer.GenomeViewer.plot`.
        Parameters
        ----------
        positions : list of int
            Genomic coordinates at which to draw marks.
        color : color_like, optional
            Color for the markers and stems. Default is "red".
        size : float, optional
            Marker size (points). Default is 8.
        line_width : float, optional
            Width of the stem lines. Default is 1.0.
        stem_length : float, optional
            Length of the stem expressed as a fraction of the axes height.
            Default is 0.06. When set to 0, the marker tip sits on the spine
            with no stem (the marker is drawn entirely above the line).
        labels : list of str, optional
            Text labels for each mark. Must be the same length as ``positions`` when
            provided. Each label is centered horizontally on its mark.
        label_rotation : float, optional
            Rotation angle for labels in degrees. Default is 90 (vertical).
        label_fontsize : float or None, optional
            Font size for labels. Defaults to the current rcParams font size.
        label_color : color_like or None, optional
            Color for labels. Defaults to ``color`` when not set.
        label_offset : float, optional
            Gap between the top of the marker and the label base, expressed as a
            fraction of the axes height. Default is 0.04.
        """
        if self._axs is None:
            raise RuntimeError(
                "You must call the `plot` method before adding axis marks."
            )
        # Filter to positions in range, keeping matching labels in sync.
        labels = list(labels)
        has_labels = len(labels) > 0
        if has_labels and len(labels) != len(positions):
            raise ValueError(
                "labels must have the same length as positions."
            )
        filtered = [
            (p, labels[i] if has_labels else None)
            for i, p in enumerate(positions)
            if self._plot_start <= p <= self._plot_end
        ]
        if not filtered:
            warn(
                "None of the provided positions are within the plotted range.",
                RuntimeWarning,
            )
            return
        ax = self._axs[0]
        effective_label_color = label_color if label_color is not None else color

        # The top spine is pushed 10 pts outward (Track._draw_track). Add the
        # same offset so y=1.0 in our mark coordinates lands on the spine line.
        from matplotlib.transforms import ScaledTranslation
        spine_offset = ScaledTranslation(
            0, 10 / 72, ax.figure.dpi_scale_trans
        )
        # x in data coords, y in axes-fraction [0=bottom, 1=top of axes box]
        transform = ax.get_xaxis_transform() + spine_offset

        spine_y = 1.0
        if stem_length > 0:
            marker_y = spine_y + stem_length
        else:
            # Scatter markers are centered on their y position. Shift upward by
            # half the marker height so the "v" tip touches the spine.
            pos = ax.get_position()
            ax_height_pts = pos.height * ax.figure.get_figheight() * ax.figure.dpi
            marker_half_frac = (size / 2.0) / ax_height_pts if ax_height_pts else 0.0
            marker_y = spine_y + marker_half_frac

        filtered_positions = [p for p, _ in filtered]
        ax.scatter(
            filtered_positions,
            [marker_y] * len(filtered_positions),
            marker="v",
            color=color,
            s=size ** 2,
            zorder=200,
            clip_on=False,
            transform=transform,
        )
        for pos, label in filtered:
            if stem_length > 0:
                ax.plot(
                    [pos, pos],
                    [spine_y, marker_y],
                    color=color,
                    linewidth=line_width,
                    transform=transform,
                    clip_on=False,
                    zorder=199,
                )
            if label is not None:
                ax.text(
                    pos,
                    marker_y + label_offset,
                    label,
                    transform=transform,
                    ha="center",
                    va="bottom",
                    rotation=label_rotation,
                    fontsize=label_fontsize,
                    color=effective_label_color,
                    clip_on=False,
                    zorder=201,
                )

        self._reserve_top_margin_for_axis_marks(
            ax=ax,
            marker_y=marker_y,
            label_offset=label_offset,
            labels=[label for _, label in filtered if label],
            label_fontsize=label_fontsize,
            label_rotation=label_rotation,
            size=size,
        )

    @staticmethod
    def _reserve_top_margin_for_axis_marks(
        ax,
        marker_y,
        label_offset,
        labels,
        label_fontsize,
        label_rotation,
        size,
    ):
        """
        Leave room above the first track for outward spines, marks, and labels.

        Without this, artists drawn above the axes bbox can be clipped until the
        figure is resized and matplotlib recomputes layout.
        """
        SPINE_OUTWARD_PT = 10
        fig = ax.figure
        fig.canvas.draw()

        pos = ax.get_position()
        fig_height_pts = fig.get_figheight() * fig.dpi
        if pos.height <= 0 or fig_height_pts <= 0:
            return

        ax_height_pts = pos.height * fig_height_pts
        extra_axes_frac = max(marker_y - 1.0, 0.0)
        extra_axes_frac += SPINE_OUTWARD_PT / ax_height_pts
        extra_axes_frac += (size / 2.0) / ax_height_pts

        if labels:
            extra_axes_frac += label_offset
            fontsize = label_fontsize or mpl.rcParams["font.size"]
            max_len = max(len(label) for label in labels)
            if label_rotation in (0, 180):
                label_pts = fontsize * max(1.2, max_len * 0.6)
            else:
                label_pts = fontsize * max(1.2, max_len * 0.55)
            extra_axes_frac += label_pts / ax_height_pts

        if ax.get_title():
            title_fontsize = ax.title.get_fontsize() or mpl.rcParams["font.size"]
            extra_axes_frac += title_fontsize * 1.4 / ax_height_pts

        extra_fig_frac = extra_axes_frac * pos.height + 0.015
        available = 1.0 - pos.y1
        if extra_fig_frac > available:
            fig.subplots_adjust(top=fig.subplotpars.top - (extra_fig_frac - available))

        fig.canvas.draw_idle()

    def _plot_group_label(self, fig, axes):
        """
        Adds a vertical group label and connecting line alongside a group of vertically stacked subplots.

        Parameters:
            fig : matplotlib.figure.Figure
                The figure object containing the subplots.
            axes : list of matplotlib.axes.Axes
                List of subplot axes, ordered top to bottom.
        """
        if len(self._group_labels.label_configs) > 0:
            fig.canvas.draw()  # Needed to update positions

        for group_config in self._group_labels.label_configs:
            # Get the bounding boxes of the top and bottom axes
            start_idx = group_config.start_track_idx
            end_idx = group_config.end_track_idx
            x = group_config.x
            x_line_offset = group_config.x_line_offset

            bbox_top = axes[start_idx].get_position()
            bbox_bottom = axes[end_idx].get_position()

            # Y coordinates in figure space
            y_top = bbox_top.y1
            y_bottom = bbox_bottom.y0
            y_center = (y_top + y_bottom) / 2

            # Match track-name (ylabel) font size and family.
            label_font = axes[start_idx].yaxis.label.get_fontproperties()

            # Add vertical label
            fig.text(
                x - x_line_offset,
                y_center,
                group_config.label,
                va="center",
                ha="center",
                rotation="vertical",
                fontproperties=label_font,
            )

            # Add vertical line beside the label
            fig.lines.append(mpl.pyplot.Line2D([x + x_line_offset, x + x_line_offset],
                                        [y_bottom, y_top],
                                        transform=fig.transFigure, color="black", linewidth=1))


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
        self._axs = None

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

        # Add group labels if specified
        self._plot_group_label(fig, axs)

        self._axs = axs
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
            heights[index] = track.layout_height()
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
