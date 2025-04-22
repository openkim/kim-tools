# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime

from kim_tools import __version__

sys.path.insert(0, os.path.abspath("../../"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "kim-tools"
copyright = (
    f"2024-{datetime.today().year} "
    "ilia Nikiforov, Ellad Tadmor, Claire Waters, "
    "Daniel Karls, Matt Bierbaum, and Eric Fuemmeler"
)
author = (
    "ilia Nikiforov, Ellad Tadmor, Claire Waters, "
    "Daniel Karls, Matt Bierbaum, and Eric Fuemmeler"
)
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery",
    "sphinx_rtd_theme",
    "sphinx.ext.todo",
]

sphinx_gallery_conf = {
    "examples_dirs": "../../examples",  # path to your example scripts
    "gallery_dirs": "auto_examples",  # path to where to save gallery generated output
}

intersphinx_mapping = {
    "ase": ("https://wiki.fysik.dtu.dk/ase/", None),
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

autoclass_content = "both"
todo_include_todos = True

autodoc_inherit_docstrings = False
autodoc_default_options = {
    "member-order": "bysource",
    "inherited-members": False,
    "show-inheritance": False,
    "private-members": True,
}

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ["_static"]

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

smartquotes = False

rst_prolog = (
    "\n.. |example_url| replace:: "
    "https://github.com/openkim-hackathons/CrystalGenomeASEExample__TD_000000654321_000"
    "\n"
)
