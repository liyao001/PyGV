# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.join(os.path.split(os.path.abspath('.'))[0], "pygv"))


# -- Project information -----------------------------------------------------

project = 'PyGV'
copyright = '2021-2024, Li Yao'
author = 'Li Yao'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.autodoc",
    # 'sphinx.ext.viewcode',
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx_gallery.gen_gallery",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to doc_source directory, that match files and
# directories to ignore when looking for doc_source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['Thumbs.db', ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

autodoc_default_options = {'autosummary': True}
image_scrapers = ('matplotlib',)
sphinx_gallery_conf = {
    "examples_dirs": "../examples",   # path to your example scripts
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated output
    "image_scrapers": image_scrapers,
    "image_srcset": ["2x"],
}
plot_formats = [('png', 300), 'pdf']
html_theme_options = {
    "show_nav_level": 4,
    # "external_links": [
    #     {"name": "Track Gallery", "url": "./auto_examples"},
    # ]
}
