"""
=============================
Coverage track for a BAM file
=============================

Create a :class:`~pygv.tracks.bam_track.CoverageTrack`.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bam_track import CoverageTrack
gv = GenomeViewer()
track = CoverageTrack(
    "../examples/data/s03.chr22.bam",
    name="RNA-seq", y_label_rotation="vertical", y_label_ha="center")
gv.add_track(track)
gv.plot("chr22", 40345142, 40367239)
