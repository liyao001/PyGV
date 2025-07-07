"""
===========
BEDPE Track
===========

"""
import matplotlib.pyplot as plt
from pygv.viewer import GenomeViewer
from pygv.tracks.bed_track import BedPETrack
from pygv.tracks.bigwig_track import BigWigTrack
gv = GenomeViewer()
interaction_track = BedPETrack(
    "../examples/data/demo.bedpe.gz",
    name="Interaction", y_label_rotation="horizontal", y_label_ha="right", height=1)
gv.add_track(interaction_track)
dnase_track = BigWigTrack(
    "../examples/data/K562_DNase_ENCFF530BKH.chr1.bw",
    name="DNase", )
gv.add_track(dnase_track)
flipped_interaction_track = BedPETrack(
    "../examples/data/demo.bedpe.gz",
    name="Interaction", y_label_rotation="horizontal", y_label_ha="right", height=1,
    flip=True
)
gv.add_track(flipped_interaction_track)
gv.plot("chr1", 155117019, 155145064)
plt.tight_layout()
