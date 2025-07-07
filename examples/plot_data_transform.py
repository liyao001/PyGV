"""
================================================
Apply data transformation to the original signal
================================================


"""
import matplotlib.pyplot as plt
import numpy as np
from pygv.viewer import GenomeViewer
from pygv.tracks.bigwig_track import PairedStrandSpecificTracks
gv = GenomeViewer()

# For this track, we don't apply any transformation
notrans_track = PairedStrandSpecificTracks(
    "../examples/data/K562_GROcap_hg38_pl.chr1.bw",
    "../examples/data/K562_GROcap_hg38_mn.chr1.bw",
    draw_y_independently=True,
    name="Original", y_label_rotation="vertical", y_label_ha="center")
gv.add_track(notrans_track)
# For this track, we apply inverse hyperbolic sine transformation
# this transformation is similar to log, but it allows for negative values
asinh_track = PairedStrandSpecificTracks(
    "../examples/data/K562_GROcap_hg38_pl.chr1.bw",
    "../examples/data/K562_GROcap_hg38_mn.chr1.bw",
    draw_y_independently=True,
    name="Asinh", y_label_rotation="vertical", y_label_ha="center")
asinh_track.data_transform = "asinh"
gv.add_track(asinh_track)
# Finally, we apply a customized function (sliding sum)
sliding_sum_track = PairedStrandSpecificTracks(
    "../examples/data/K562_GROcap_hg38_pl.chr1.bw",
    "../examples/data/K562_GROcap_hg38_mn.chr1.bw",
    draw_y_independently=True,
    name="Customized\nfunc", y_label_rotation="vertical", y_label_ha="center")
# Define a function to calculate the sliding sum (window size: 49)
# the following function requires numpy version >= 1.20.0
sliding_sum_track.data_transform = lambda x: np.pad(
    np.lib.stride_tricks.sliding_window_view(x, 49).sum(axis=-1), (24, 24))
gv.add_track(sliding_sum_track)
gv.plot("chr1", 201954851, 201955948)
plt.tight_layout()
