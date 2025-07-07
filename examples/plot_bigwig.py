"""
============
BigWig Track
============

Create a simple :class:`~pygv.tracks.bigwig_track.BigWigTrack` Track.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bigwig_track import BigWigTrack
gv = GenomeViewer()
track = BigWigTrack(
    "../examples/data/K562_DNase_hg38_ENCFF413AHU.chr19.bigWig",
    name="DNase", )
gv.add_track(track)
bar_track = BigWigTrack(
    "../examples/data/K562_DNase_hg38_ENCFF413AHU.chr19.bigWig",
    name="DNase", plot_type="bar")
gv.add_track(bar_track)
gv.plot("chr19", 7496780, 7498169)
