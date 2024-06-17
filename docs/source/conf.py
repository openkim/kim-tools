# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'kim-test-utils'
copyright = '2024, ilia Nikiforov, Ellad Tadmor, and Eric Fuemmeler'
author = 'ilia Nikiforov, Ellad Tadmor, and Eric Fuemmeler'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc","sphinx.ext.napoleon","sphinx_autodoc_typehints","sphinx.ext.intersphinx","sphinx.ext.todo"]

intersphinx_mapping = {'ase': ('https://wiki.fysik.dtu.dk/ase/', None)}

autoclass_content="both"
todo_include_todos = True

autodoc_inherit_docstrings = False
autodoc_default_options = {
    'member-order': 'bysource',
    'inherited-members': False,
    'show-inheritance': False
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
