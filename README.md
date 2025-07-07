# Python Genome Viewer (PyGV)

PyGV, short for *Py*thon *G*enome *V*iewer, is a powerful and user-friendly Python package designed to empower
researchers and genomics enthusiasts with the tools they need to create publication-ready, highly customizable genomic
visualizations seamlessly integrated into their daily workflows.

## Key Features and Advantages

* **Everything in Python**: PyGV allows you to effortlessly configure genomic tracks directly within Python code.
  Each track is represented as a Python class, making it easy to fine-tune track properties, layouts, and styles to suit
  your specific research needs.
  As PyGV is built on top of Matplotlib, you can further harness the familiar syntax used for Matplotlib,
  giving you the power to customize your genomic visualizations to the fullest extent.
* **Flexible data filtering**: Users can define what to plot by applying customized filters during signal reading.
  For example, when plotting alignments, users can use Python's lambda expressions to selectively display specific
  reads, such as the first read or reads mapped to the forward strand. 
  This level of customization extends to genome annotations as well, allowing users to effortlessly select which 
  features (e.g. transcripts) to display with just one line of code. This capability is useful when dealing 
  with extensive genomic data, enabling researchers to focus on the specific elements relevant to their analyses 
  without being overwhelmed by unnecessary details.
* **Rich Genomic Data Support**: PyGV offers support for a wide range of genomic data formats, allowing you to
  visualize genes, regulatory elements, sequencing alignments, genomic interactions, and other genomics features in the
  way that best serves your research. See a gallery of all supported tracks [here]().
* **High-resolution Figures**: PyGV is committed to producing high-quality figures. Our goal is to minimize the
  additional effort you might otherwise spend adjusting fonts, colors, and other details in software like Adobe
  Illustrator. With PyGV, you can create polished, high-quality figures directly from your genomic visualizations,
  saving you time and ensuring that your visualizations are ready to meet the different standards of publication.
* **Reproducibility and Collaboration**: With PyGV, creating and sharing reproducible genomic visualizations becomes
  a straightforward task, essential for research collaboration and documentation.

## Get Started

Install PyGV with pip:

```shell
pip install pygv
```

Now let's have some fun:

```python
from pygv.viewer import GenomeViewer
from pygv.tracks import gtf_track

gv = GenomeViewer()
gencode_track = gtf_track.GtfTrack(
    "~/gencode.v34lift37.annotation.sorted.gtf.gz",
    name="GENCODE", show_genes=True, show_transcript_id=True,
    filters=lambda x: x.transcript_id in {"ENST00000332995.11_1", "ENSG00000112137.17_4",
                                          "ENST00000379350.5_1", "ENST00000379335.7_1"},
    annotation_formatter=lambda x: x.split(".")[0]
)
gv.add_track(gencode_track)
gv.plot("chr6", 12714999, 13292716)
plt.show()
```