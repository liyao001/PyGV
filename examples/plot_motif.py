"""
===========
Motif Track
===========

:class:`~pygv.tracks.motif_track.MotifTrack` searches a reference FASTA for a
nucleotide motif and draws each hit as a stranded interval, similar to IGV's
`Motif Finder <https://igv.org/doc/desktop/#UserGuide/tools/motif_finder/>`_.

Patterns can be a literal sequence, IUPAC ambiguity codes (for example
``AAAAAARNR``), or a nucleotide regular expression (for example
``TATA[AT]A[AT]A`` or ``TATAAA+``).
"""
import matplotlib.pyplot as plt
from pygv.viewer import GenomeViewer
from pygv.tracks.motif_track import MotifTrack

gv = GenomeViewer()

# %%
# :meth:`~pygv.tracks.motif_track.MotifTrack.pair` adds the IGV-style pair of
# tracks: plus-strand hits in blue and minus-strand hits in red. Overlapping
# matches stay on separate lanes because the default show mode is expanded.
plus, minus = MotifTrack.pair(
    "../examples/data/motif_demo.fa",
    "TATA[AT]A[AT]A",
    name="TATA box",
)
gv.add_tracks((plus, minus))

# %%
# IUPAC codes are expanded the same way as in IGV. This pattern is six adenines,
# a purine, any base, and another purine.
iupac = MotifTrack(
    "../examples/data/motif_demo.fa",
    "AAAAAARNR",
    strand="+",
    name="AAAAAARNR",
)
gv.add_track(iupac)

# %%
gv.plot("chr1", 0, 250, height_scale_factor=0.6)
plt.tight_layout()
