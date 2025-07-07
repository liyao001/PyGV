"""
===================
Dynamic value track
===================

In this example, we will use PyGV to plot the sine function in a certain range.
"""
import matplotlib.pyplot as plt
import numpy as np
from pygv.viewer import GenomeViewer
from pygv.tracks.track import DynamicValueTrack
gv = GenomeViewer()
track = DynamicValueTrack("", name="Demo")
gv.add_track(track)

demo_start = 10000
demo_end = 10050
track.values = np.sin(np.arange(demo_end-demo_start))
gv.plot("chr19", demo_start, demo_end)
plt.show()
