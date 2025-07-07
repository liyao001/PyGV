"""
=============================================================
Read alignment track with split parts joined by thinner lines
=============================================================

Create a :class:`~pygv.tracks.bam_track.SplicedReadTrack` to show the
alignment in a BAM file, the spliced parts will be rendered as thinner lines.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bam_track import SplicedReadTrack
gv = GenomeViewer()
track = SplicedReadTrack(
    "../examples/data/s03.chr22.bam",
    name="RNA-seq")
track.sampling_ratio = 0.3
track.color_reads_by = "first of pair strand"
track.line_width = 0.3
gv.add_track(track)
gv.plot("chr22", 40349785, 40353365)
