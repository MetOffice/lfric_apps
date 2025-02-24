# Configuration file for the Sphinx documentation builder.
#
# -- imports -----------------------------------------------------------------
import os

import docutils


simulations-systems

# -- Project information -----------------------------------------------------

project = 'LFRic Apps'
author = 'Simulation IT'
copyright = '2025 Met Office. All rights reserved'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_sitemap',
    'sphinx_design'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

html_static_path = ["_static"]
html_css_files = ["custom.css"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'
# html_title = "LFRic Apps"

# Generate the sitemap info, this will need updating when we have versioned docs
html_baseurl = os.environ.get("SPHINX_HTML_BASE_URL", "https://metoffice.github.io/lfric_apps")
sitemap_locales = [None]
sitemap_url_scheme = "{link}"

# Hide the link which shows the rst markup
html_show_sourcelink = False

html_theme_options = {
    "navigation_with_keys": True,
    "use_edit_page_button": True,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_align": "content",
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/MetOffice/lfric_apps",
            "icon": "fa-brands fa-github"
        },
        {
            "name": "GitHub Discussions",
            "url": "https://github.com/MetOffice/simulation-systems/discussions",
            "icon": "far fa-comments",
        }
    ],
    "logo": {
        "text": "LFRic Apps",
        "image_light": "_static/MO_SQUARE_black_mono_for_light_backg_RBG.png",
        "image_dark": "_static/MO_SQUARE_for_dark_backg_RBG.png",
    },
    "secondary_sidebar_items": {
        "**/*": ["page-toc", "edit-this-page", "show-glossary"],
        "index": [],
    },
    "footer_start": ["crown-copyright"],
    "footer_center": ["sphinx-version"],
    "footer_end": ["theme-version"],
    "primary_sidebar_end": []
}

html_sidebars = {
    "index": []
}

# Provides the Edit on GitHub link in the generated docs.
html_context = {
    "display_github": True,
    "github_user": "MetOffice",
    "github_repo": "lfric_apps",
    "github_version": "main",
    "doc_path": "documentation/source"
}

# Enable numbered references to e.g. figures.
#
numfig = True

# Exclude files from Sphinx processing
exclude_patterns = ['common_links.rst']

# Create rst_epilog variable to allow concatenation of epilog parts which will
# be included at the end of every rst file.
rst_prolog = ""

# Add the contents of the common_links file to the epilog.
with open('common_links.rst') as file:
    rst_prolog += file.read()


def superscript_substitution_role(name, rawtext, text, lineno, inliner,
                                  options={}, content=[]):
    node = docutils.nodes.superscript()
    node2 = docutils.nodes.substitution_reference(refname=text)
    node += [node2]
    return [node], []


def setup(app):
    app.add_role('superscript_substitution', superscript_substitution_role)
