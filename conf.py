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
sys.path.insert(0, os.path.abspath(''))


# -- Project information -----------------------------------------------------

project = 'Process Engineering Project'
copyright = '2021, Antonio Rocha Azevedo'
author = 'Antonio Rocha Azevedo'

# The full version, including alpha/beta/rc tags
# release = '2021'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon',
              'sphinx.ext.autosectionlabel', 'nbsphinx',
              'sphinx.ext.viewcode', 'nbsphinx_link'
]

autodoc_member_order = 'bysource'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Highlighting -----------------------------------------------------------

from pygments.styles import get_all_styles
from pygments.styles import get_style_by_name

# styles = list(get_all_styles())
# print(styles)

pygments_style = 'rainbow_dash'

# -- Options for HTML output -------------------------------------------------

# html_logo = '_assets/logo.gif'
html_favicon = '_assets/favicon.ico'

html_theme = 'furo'
html_theme_options = {
    "light_css_variables": {
        # "color-brand-primary": "red",                 # Links in the side-bar
        # "color-brand-content": "#CC3333",             # function names etc.
        # "color-admonition-background": "00FF00",
    },
    "dark_css_variables": {
        "color-brand-primary": "#C9FFEE",
        "color-brand-content": "#00FF6A",
        # "color-admonition-background": "00FF00",
    },
    # "announcement": "Welcome to ProcessEngineering's documentation!"
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']