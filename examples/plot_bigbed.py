"""
============
BigBed Track
============

Create a :class:`~pygv.tracks.bigbed_track.BigBed6Track`.
"""
import matplotlib.pyplot as plt
from pygv.viewer import GenomeViewer
from pygv.tracks.bigbed_track import BigBed6Track
gv = GenomeViewer()
track = BigBed6Track(
    "../examples/data/test.bigBed",
    name="With names", )
gv.add_track(track)
no_name_track = BigBed6Track(
    "../examples/data/test.bigBed",
    name="No names", )
no_name_track.show_name = False
gv.add_track(no_name_track)
gv.plot("chr14", 100235930, 100242295)
plt.tight_layout()
