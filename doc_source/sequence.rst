Tracks for sequences
********************

Sequence Logo
-------------
.. automodule:: pygv.tracks.logo_track
    :members:

Motif hits
----------
Search a reference FASTA for a nucleotide motif and draw the hits as stranded
interval features, similar to IGV's
`Motif Finder <https://igv.org/doc/desktop/#UserGuide/tools/motif_finder/>`_.

Patterns may be a literal sequence (``ACCGCT``), IUPAC ambiguity codes
(``AAAAAARNR``), or a nucleotide regular expression (``TATA[AT]A[AT]A``,
``TATAAA+``). Use :meth:`~pygv.tracks.motif_track.MotifTrack.pair` to add the
IGV-style plus (blue) and minus (red) tracks together.

.. automodule:: pygv.tracks.motif_track
    :members:
