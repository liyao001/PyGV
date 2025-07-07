"""
===================
Sequence Logo Track
===================

A sequence logo is a graphical representation of the
sequence conservation and variation in a set of aligned sequences,
typically representing a group of related DNA, RNA, or protein sequences.
The primary purpose of a sequence logo is to visualize the sequence motifs
or patterns that are conserved across the aligned sequences.

In PyGV, you can create a :class:`~pygv.tracks.logo_track.LogoTrack` to show Sequence Logos.
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pygv.viewer import GenomeViewer
from pygv.tracks.track import DynamicValueTrack
from pygv.tracks.logo_track import LogoTrack


gv = GenomeViewer()

track = LogoTrack("", name="logo")
gv.add_track(track)

dyn_track = DynamicValueTrack("", name="data")
gv.add_track(dyn_track)

demo_start = 10000
demo_end = 10129
track.values = pd.read_csv("../examples/data/saliency.csv", index_col=0)
dyn_track.values = np.sin(np.arange(demo_end-demo_start))
gv.plot("chr19", demo_start, demo_end)
plt.tight_layout()
plt.show()
