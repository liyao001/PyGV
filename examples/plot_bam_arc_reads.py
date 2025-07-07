"""
====================================================
Read alignment track with split parts joined by arcs
====================================================

Create a :class:`~pygv.tracks.bam_track.ReadArcTrack` to show spliced reads from a BAM file.
In contrast to `SplicedReadTrack`, this track shows the spliced parts as arcs.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bam_track import ReadArcTrack
gv = GenomeViewer()
track = ReadArcTrack(
    "../examples/data/s03.chr22.bam",
    name="RNA-seq", )
track.sampling_ratio = 0.3
gv.add_track(track)
gv.plot("chr22", 40349785, 40353365)
