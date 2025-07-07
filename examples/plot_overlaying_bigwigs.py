"""
========================
Overlaying BigWig Tracks
========================

Create a :class:`~pygv.tracks.bigwig_track.OverlayingTrack`
to show signals from multiple BigWig files within a single track.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bigwig_track import OverlayingTrack
gv = GenomeViewer()
track = OverlayingTrack(
    ("../examples/data/K562_H3K27ac_ENCFF779QTH.chr19.bigWig",
     "../examples/data/K562_DNase_hg38_ENCFF413AHU.chr19.bigWig",),
    labels=("H3K27ac", "DNase", ))
gv.add_track(track)
gv.plot("chr19", 13232024, 13237758)
