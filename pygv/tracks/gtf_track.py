import math
import os
import re
from collections import OrderedDict, namedtuple
from typing import Any, Optional

import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from pydantic import Field, PrivateAttr

from .bed_track import _LaneRegistry
from .track import AnnotationTrack
from .types import FilterFn


class GtfTrack(AnnotationTrack):
    """GTF annotation track."""

    filters: FilterFn = Field(
        default=None,
        description="Callable returning True/False to keep/discard GTF records",
    )
    show_genes: bool = Field(
        default=False,
        description="Whether to display gene-level features in addition to transcripts",
    )
    annotation_formatter: Any = Field(
        default=None,
        description="Callable used to format gene/transcript labels",
    )
    show_transcript_id: bool = Field(
        default=False,
        kw_only=True,
        description="Show transcript IDs instead of gene names when available",
    )

    _gtf_file: Optional[str] = PrivateAttr(default=None)
    _fields: tuple = PrivateAttr(default=())
    _gtf_obj: Any = PrivateAttr(default=None)
    _GtfRecord: Any = PrivateAttr(default=None)
    _parser: Any = PrivateAttr(default=None)

    def __init__(
        self,
        track: str,
        filters: FilterFn = None,
        show_genes: bool = False,
        annotation_formatter: Any = None,
        **data: Any,
    ) -> None:
        super().__init__(
            track=track,
            filters=filters,
            show_genes=show_genes,
            annotation_formatter=annotation_formatter,
            **data,
        )

    def _get(self, chromosome, start, end):
        if self._parser is None:
            return
        yield from self._parser(chromosome, start, end)

    def _pysam_parser(self, chromosome, start, end):
        import pysam

        raw_hits = dict()
        regex = r"(.*?)\s\"(.*?)\";"
        try:
            for i, hit in enumerate(
                self._gtf_obj.fetch(chromosome, start, end, parser=pysam.asGTF())
            ):
                if hit.feature == "transcript" or (
                    self.show_genes and hit.feature == "gene"
                ):
                    if self.filters is not None:
                        if hit.feature == "gene":
                            hit.transcript_id = hit.gene_id
                        if not self.filters(hit):
                            continue
                    if hit.feature == "gene":
                        k = f"{hit.gene_id}|0"
                    else:
                        k = f"{hit.gene_id}|{hit.transcript_id}"
                    attributes = {}
                    matches = re.finditer(regex, hit.attributes, re.MULTILINE)
                    for matchNum, match in enumerate(matches, start=1):
                        key, value = match.groups()
                        attributes[key.strip()] = value
                    if k not in raw_hits:
                        raw_hits[k] = {
                            "contig": hit.contig,
                            "feature": hit.feature,
                            "doc_source": hit.source,
                            "start": hit.start,
                            "end": hit.end,
                            "score": hit.score,
                            "strand": hit.strand,
                            "frame": hit.frame,
                            "attributes": attributes,
                            "transcript_id": hit.transcript_id,
                            "gene_id": hit.gene_id,
                            "exons": [],
                        }
                    else:
                        raw_hits[k]["contig"] = hit.contig
                        raw_hits[k]["feature"] = hit.feature
                        raw_hits[k]["doc_source"] = hit.source
                        raw_hits[k]["start"] = hit.start
                        raw_hits[k]["end"] = hit.end
                        raw_hits[k]["score"] = hit.score
                        raw_hits[k]["strand"] = hit.strand
                        raw_hits[k]["frame"] = hit.frame
                        raw_hits[k]["attributes"] = attributes
                        raw_hits[k]["transcript_id"] = hit.transcript_id
                        raw_hits[k]["gene_id"] = (hit.gene_id,)
                elif hit.feature == "exon":
                    if self.filters is not None:
                        if not self.filters(hit):
                            continue
                    k = f"{hit.gene_id}|{hit.transcript_id}"
                    if k in raw_hits:
                        raw_hits[k]["exons"].append((hit.start, hit.end))
                    else:
                        raw_hits[k] = {"exons": [(hit.start, hit.end)]}
            if self.show_genes:
                # order by: gene ID, gene or transcript, segment length, start position
                sorted_keys = sorted(
                    raw_hits,
                    key=lambda k: (
                        raw_hits[k]["gene_id"][0],
                        raw_hits[k]["transcript_id"] == raw_hits[k]["gene_id"],
                        1 / (raw_hits[k]["end"] - raw_hits[k]["start"]),
                        raw_hits[k]["start"],
                    ),
                )
                hits = OrderedDict({k: raw_hits[k] for k in sorted_keys})
                raw_hits = hits
            for transcript_id, transcript in raw_hits.items():
                try:
                    yield self._GtfRecord(**transcript)
                except Exception as e:
                    print(e)
                    pass
        except ValueError:
            # in case no feature is available in that window
            return

    def _pd_parser(self, chromosome, start, end):
        df = self._gtf_obj
        if df is None or len(df) == 0:
            return
        window = df[
            (df["seqname"] == chromosome)
            & (df["start"] <= end)
            & (df["end"] > start)
        ]
        raw_hits = {}
        regex = r"(.*?)\s\"(.*?)\";"
        for row in window.itertuples(index=False):
            attributes = {}
            for match in re.finditer(regex, str(row.attribute), re.MULTILINE):
                key, value = match.groups()
                attributes[key.strip()] = value
            gene_id = attributes.get("gene_id", "")
            transcript_id = attributes.get("transcript_id", "")
            start0 = int(row.start) - 1
            end0 = int(row.end)
            if row.feature == "transcript" or (
                self.show_genes and row.feature == "gene"
            ):
                if self.filters is not None and not self.filters(row):
                    continue
                k = (
                    f"{gene_id}|0"
                    if row.feature == "gene"
                    else f"{gene_id}|{transcript_id}"
                )
                exons = list(raw_hits.get(k, {}).get("exons", []))
                raw_hits[k] = {
                    "contig": row.seqname,
                    "feature": row.feature,
                    "doc_source": row.source,
                    "start": start0,
                    "end": end0,
                    "score": row.score,
                    "strand": row.strand,
                    "frame": row.frame,
                    "attributes": attributes,
                    "transcript_id": transcript_id,
                    "gene_id": gene_id,
                    "exons": exons,
                }
            elif row.feature == "exon":
                if self.filters is not None and not self.filters(row):
                    continue
                k = f"{gene_id}|{transcript_id}"
                if k in raw_hits:
                    raw_hits[k]["exons"].append((start0, end0))
                else:
                    raw_hits[k] = {"exons": [(start0, end0)]}
        for transcript in raw_hits.values():
            if "contig" not in transcript:
                continue
            try:
                yield self._GtfRecord(**transcript)
            except Exception as e:
                print(e)

    def model_post_init(self, __context: Any) -> None:
        super().model_post_init(__context)
        if not os.path.exists(self.track) and not self.track.startswith("http"):
            raise ValueError

        self._gtf_file = self.track
        self._fields = (
            "contig",
            "feature",
            "doc_source",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "gene_id",
            "transcript_id",
            "attributes",
            "exons",
        )
        use_pysam = 1
        try:
            import pysam
        except ImportError:
            use_pysam = 0

        if use_pysam and self.track.endswith(".gtf.gz"):
            if os.path.exists(self.track + ".tbi"):
                use_pysam = 1
            elif self.track.startswith("http"):
                use_pysam = 1
            else:
                try:
                    pysam.tabix_index(self.track)
                    use_pysam = 1
                except Exception as e:
                    use_pysam = 0
                    print("Failed to index gtf file", self.track, e)
        else:
            use_pysam = 0
        if use_pysam:
            self._gtf_obj = pysam.TabixFile(self.track)
            self._parser = self._pysam_parser
        else:
            import pandas as pd

            self._gtf_obj = pd.read_csv(
                self.track,
                sep="\t",
                header=None,
                comment="#",
                names=[
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "attribute",
                ],
            )
            self._parser = self._pd_parser

        self._GtfRecord = namedtuple("GtfRecord", self._fields)
        self._plot_block = 1
        if self.annotation_formatter is None or not callable(self.annotation_formatter):
            self.annotation_formatter = lambda x: x

    def _pre_plot_hook(self, chromosome, start, end, **kwargs):
        """
        Build non-overlapping tracks

        Parameters
        ----------
        chromosome : str
            Chromosome
        start : int
            start of visible window
        end : int
            end of visible window
        Returns
        -------

        """
        # clean up lane registry
        self._lane_registries = []
        added = set()
        for interval in self._get(chromosome=chromosome, start=start, end=end):
            active_lane = None
            start_loc = int(interval.start)
            end_loc = int(interval.end)
            visible_start = max(start_loc, start)
            visible_end = min(end_loc, end)

            if self.hide_visual_dup:
                k = (visible_start, visible_end, interval.strand)
                if k in added:
                    continue
                added.add(k)

            if len(self._lane_registries) == 0:
                self._lane_registries.append(_LaneRegistry())

            for lr in self._lane_registries:
                if lr.max_coord is None:
                    active_lane = lr.offset
                    lr.min_coord = start_loc
                    lr.max_coord = end_loc
                    lr.features.append(interval)
                    break
                else:
                    if lr.max_coord < start_loc - self.padding_left:
                        active_lane = lr.offset
                        lr.min_coord = min(lr.min_coord, start_loc)
                        lr.max_coord = max(lr.max_coord, end_loc)
                        lr.features.append(interval)
                        break

            if (
                type(self.allowed_feature_lanes) is int
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

    def _draw_track(self, chromosome, start, end, ax, index=1, **kwargs):
        super(GtfTrack, self)._draw_track(
            chromosome=chromosome, start=start, end=end, ax=ax, index=index, **kwargs
        )
        self._ax.set_xlim((start, end))

        self._small_relative = 0.004 * (end - start)
        for lane in self._lane_registries:
            for interval in lane.features:
                start_loc = int(interval.start)
                end_loc = int(interval.end)
                visible_start = max(start_loc, start)
                visible_end = min(end_loc, end)
                active_lane = lane.offset

                real_active_line = (self.patch_height + self.lane_space) * active_lane

                if not self._plot_block:
                    rec = Rectangle(
                        xy=(
                            start_loc,
                            -1 * real_active_line - (self.patch_height / 2),
                        ),
                        width=end_loc - start_loc,
                        height=self.patch_height,
                        facecolor=self.color,
                        linewidth=self.line_width,
                    )
                    self._ax.add_patch(rec)
                else:
                    # init
                    exon_starts = list(
                        map(int, [coords[0] for coords in interval.exons])
                    )
                    exon_ends = list(map(int, [coords[1] for coords in interval.exons]))
                    if len(exon_starts) == len(exon_ends):
                        if len(exon_starts) == 0:
                            exon_starts.append(start_loc)
                            exon_ends.append(end_loc)
                        self._ax.plot(
                            (
                                start_loc if start_loc >= start else visible_start,
                                end_loc - 1 if end_loc <= end else visible_end,
                            ),
                            (-1 * real_active_line, -1 * real_active_line),
                            color=self.line_color,
                            linewidth=self.line_width,
                            clip_on=True,
                        )

                        # plot small arrows over the backbone
                        if end_loc - start_loc > self._small_relative:
                            pos = np.arange(
                                visible_start + self._small_relative,
                                visible_end + self._small_relative,
                                int(self.arrow_interval * self._small_relative),
                            )
                            for xpos in pos:
                                self._plot_gene_direction(
                                    ax, xpos, -1 * real_active_line, interval.strand
                                )

                        patches = []
                        for i, (s, e) in enumerate(zip(exon_starts, exon_ends)):
                            if s < end or e > start:
                                # in case overflow
                                adjusted_size = e - s
                                plot_end = e
                                if plot_end > visible_end:
                                    adjusted_size -= plot_end - visible_end
                                p = Rectangle(
                                    xy=(
                                        s,
                                        -1 * real_active_line
                                        - (self.patch_height / 2),
                                    ),
                                    width=adjusted_size,
                                    clip_on=True,
                                    height=self.patch_height,
                                )
                                patches.append(p)

                        self._ax.add_collection(
                            PatchCollection(
                                patches,
                                edgecolors=self.edge_color,
                                facecolors=self.color,
                                linewidths=self.line_width,
                                zorder=100,
                                clip_on=True,
                            )
                        )
                    else:
                        rec = Rectangle(
                            xy=(
                                start_loc,
                                -1 * real_active_line - self.patch_height / 2,
                            ),
                            width=end_loc - start_loc,
                            height=self.patch_height,
                            edgecolor=self.edge_color,
                            facecolor=self.color,
                            linewidth=self.line_width,
                            **kwargs,
                        )
                        self._ax.add_patch(rec)

                cond_a = "gene_name" in interval.attributes and self.show_name
                cond_b = (
                    "transcript_id" in interval.attributes and self.show_transcript_id
                )
                if cond_a or cond_b:
                    cond_b = False if interval.feature == "gene" else cond_b
                    if start_loc > start and interval.strand == "+":
                        self._ax.text(
                            x=start_loc - self._small_relative,
                            y=-1 * real_active_line,
                            color=self.font_color,
                            size=self.font_size,
                            s=self.annotation_formatter(
                                interval.attributes["transcript_id"]
                                if cond_b
                                else interval.attributes["gene_name"]
                            ),
                            ha="right",
                            va="center",
                            clip_on=False,
                            zorder=101,
                        )
                    elif end_loc < end and interval.strand == "-":
                        self._ax.text(
                            x=end_loc + self._small_relative,
                            y=-1 * real_active_line,
                            color=self.font_color,
                            size=self.font_size,
                            s=self.annotation_formatter(
                                interval.attributes["transcript_id"]
                                if cond_b
                                else interval.attributes["gene_name"]
                            ),
                            ha="left",
                            va="center",
                            clip_on=False,
                            zorder=101,
                        )
                    else:
                        self._ax.text(
                            x=(visible_end + visible_start) / 2,
                            y=-1 * real_active_line,
                            color=self.font_color,
                            size=self.font_size,
                            s=self.annotation_formatter(
                                interval.attributes["transcript_id"]
                                if cond_b
                                else interval.attributes["gene_name"]
                            ),
                            ha="center",
                            va="center",
                            clip_on=False,
                            bbox=dict(
                                boxstyle="round",
                                fc="w",
                                alpha=self.font_box_alpha,
                                lw=0.1,
                            ),
                            zorder=101,
                        )

        self._ax.set_yticks([])
        # remove minor ticks
        self._ax.set_yticks([], minor=True)

        if self._plot_block:
            n = len(self._lane_registries)
            units = math.ceil(n / self.features_per_lane)
            units = units if units >= 1 else 1
            ylim_lower = -(self.patch_height + self.lane_space) * len(
                self._lane_registries
            )
            aesthetic_lower = (
                -(self.patch_height + self.lane_space)
                * self.features_per_lane
                * units
            )
            ylim_lower = ylim_lower if ylim_lower < aesthetic_lower else aesthetic_lower

            self._ax.set_ylim((ylim_lower, self.patch_height / 2 + 0.05))
        else:
            n = len(self._lane_registries)
            self._ax.set_ylim((-1.5, 1.5))

        if index != 0:
            # remove major ticks
            self._ax.set_xticks([])
            # remove minor ticks
            self._ax.set_xticks([], minor=True)
            # self.ax.margins(0)
