"""
============
BigWig Track
============

Create a simple BigWig Track.
"""
from pygv.viewer import GenomeViewer
from pygv.tracks.bigwig_track import BigWigTrack
gv = GenomeViewer()
track = BigWigTrack(
    "../examples/data/K562_DNase_ENCFF530BKH.chr1.bw",
    name="DNase", )
gv.add_track(track)
gv.plot("chr1", 73917, 78917)
