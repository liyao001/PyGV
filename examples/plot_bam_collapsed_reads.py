"""
=====================
Collapsed reads track
=====================

Create a :class:`~pygv.tracks.bam_track.CollapsedReadTrack` to show aligned reads in a BAM file.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bam_track import CollapsedReadTrack
gv = GenomeViewer()
track = CollapsedReadTrack(
    "../examples/data/s03.chr22.bam",
    name="RNA-seq", )
track.sampling_ratio = 0.3
gv.add_track(track)
gv.plot("chr22", 40345142, 40367239)
