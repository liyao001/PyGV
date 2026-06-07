"""
==========
Axis Marks
==========

Add triangle marks with optional text labels to the top coordinate
spine using :meth:`~pygv.viewer.GenomeViewer.set_axis_marks`.
"""
import matplotlib.pyplot as plt
from pygv.viewer import GenomeViewer
from pygv.tracks.bigwig_track import BigWigTrack

gv = GenomeViewer()
dnase_track = BigWigTrack(
    "../examples/data/K562_DNase_ENCFF530BKH.chr1.bw",
    name="DNase",
)
gv.add_track(dnase_track)
gv.plot("chr1", 155117019, 155145064)

# Add marks at specific genomic positions with centered text labels.
gv.set_axis_marks(
    positions=[155121528, 155127872, 155135811, 155140407],
    labels=["rs001", "rs002", "rs003", "rs004"],
    color="red",
    size=8,
    line_width=1.0,
    stem_length=0.0,
    label_rotation=0,
    label_offset=0.04,
)
