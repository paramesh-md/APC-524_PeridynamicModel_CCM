# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os

sys.path.insert(0, os.path.abspath('.'))

project = 'PyPeridynamics'
copyright = '2023, Parameshwaran Pasupathy'
author = 'Parameshwaran Pasupathy'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser",
              "sphinx.ext.autodoc",
              "sphinx.ext.napoleon",
              "sphinx.ext.viewcode",
              "sphinx.ext.intersphinx",
              "sphinx.ext.mathjax",
              "sphinx_autodoc_typehints",
              "sphinx_copybutton",]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.venv', '.nox', '.pytest_cache', 'lib']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

myst_enable_extensions = [
    "colon_fence",
    ]
