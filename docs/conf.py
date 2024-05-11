"""
Configuration file for the Sphinx documentation builder.
"""

import os

# This needs to be done before aiapy or sunpy is imported
os.environ["PARFIVE_HIDE_PROGRESS"] = "True"

import datetime  # NOQA: E402
import warnings  # NOQA: E402
from pathlib import Path  # NOQA: E402

from astropy.utils.exceptions import AstropyDeprecationWarning  # NOQA: E402
from matplotlib import MatplotlibDeprecationWarning  # NOQA: E402
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyPendingDeprecationWarning  # NOQA: E402
from sunpy_sphinx_theme import PNG_ICON  # NOQA: E402

from aiapy import __version__  # NOQA: E402

# -- Project information -------------------------------------------------------
project = "aiapy"
author = "AIA Instrument Team"
copyright = f"{datetime.datetime.now(datetime.timezone.utc).year}, {author}"  # NOQA: A001
release = __version__
is_development = ".dev" in __version__

# Need to make sure that our documentation does not raise any of these
warnings.filterwarnings("error", category=SunpyDeprecationWarning)
warnings.filterwarnings("error", category=SunpyPendingDeprecationWarning)
warnings.filterwarnings("error", category=MatplotlibDeprecationWarning)
warnings.filterwarnings("error", category=AstropyDeprecationWarning)

linkcheck_ignore = [
    r"https://doi.org/\d+",
    r"https://element.io/\d+",
    r"https://github.com/\d+",
    r"https://docs.sunpy.org/\d+",
]
linkcheck_anchors = False

# -- General configuration -----------------------------------------------------
extensions = [
    "matplotlib.sphinxext.plot_directive",
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_changelog",
    "sphinx_gallery.gen_gallery",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sunpy.util.sphinx.doctest",
    "sunpy.util.sphinx.generate",
    "sphinxext.opengraph",
    "sphinx_design",
    "sphinx_copybutton",
    "hoverxref.extension",
]
automodapi_toctreedirnm = "generated/api"
html_extra_path = ["robots.txt"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"
default_role = "obj"
napoleon_use_rtype = False
napoleon_google_docstring = False
napoleon_use_param = False
suppress_warnings = ["app.add_directive"]
nitpicky = True
# This is not used. See docs/nitpick-exceptions file for the actual listing.
nitpick_ignore = []
with Path("nitpick-exceptions").open() as nitpick_exceptions:
    for line in nitpick_exceptions:
        if line.strip() == "" or line.startswith("#"):
            continue
        dtype, target = line.split(None, 1)
        target = target.strip()
        nitpick_ignore.append((dtype, target))

# -- Options for sphinxext-opengraph ------------------------------------------
ogp_image = "https://raw.githubusercontent.com/sunpy/sunpy-logo/master/generated/sunpy_logo_word.png"
ogp_use_first_image = True
ogp_description_length = 160
ogp_custom_meta_tags = [
    '<meta property="og:ignore_canonical" content="true" />',
]

# -- Options for sphinx-copybutton ---------------------------------------------
# Python Repl + continuation, Bash, ipython and qtconsole + continuation, jupyter-console + continuation
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

# -- Options for hoverxref -----------------------------------------------------
if os.environ.get("READTHEDOCS"):
    hoverxref_api_host = "https://readthedocs.org"
    if os.environ.get("PROXIED_API_ENDPOINT"):
        # Use the proxied API endpoint
        # - A RTD thing to avoid a CSRF block when docs are using a
        #   custom domain
        hoverxref_api_host = "/_"
hoverxref_tooltip_maxwidth = 600  # RTD main window is 696px
hoverxref_auto_ref = True
hoverxref_mathjax = True
# hoverxref has to be applied to these
hoverxref_domains = ["py"]
hoverxref_role_types = {
    # roles with py domain
    "attr": "tooltip",
    "class": "tooltip",
    "const": "tooltip",
    "data": "tooltip",
    "exc": "tooltip",
    "func": "tooltip",
    "meth": "tooltip",
    "mod": "tooltip",
    "obj": "tooltip",
    # roles with std domain
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
}

# -- Options for intersphinx extension -----------------------------------------
intersphinx_mapping = {
    "astropy": ("https://docs.astropy.org/en/stable/", None),
    "cupy": ("https://docs.cupy.dev/en/stable/", None),
    "drms": ("https://docs.sunpy.org/projects/drms/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "numpy": ("https://numpy.org/doc/stable/", (None, "http://www.astropy.org/astropy-data/intersphinx/numpy.inv")),
    "parfive": ("https://parfive.readthedocs.io/en/stable/", None),
    "pyerfa": ("https://pyerfa.readthedocs.io/en/stable/", None),
    "python": ("https://docs.python.org/3/", (None, "http://www.astropy.org/astropy-data/intersphinx/python3.inv")),
    "reproject": ("https://reproject.readthedocs.io/en/stable/", None),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://www.astropy.org/astropy-data/intersphinx/scipy.inv"),
    ),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable/", None),
}

# -- Options for HTML output ---------------------------------------------------
html_theme = "sunpy"
graphviz_output_format = "svg"
graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]

# -- Sphinx Gallery ------------------------------------------------------------
# JSOC email os env
# see https://github.com/sunpy/sunpy/wiki/Home:-JSOC
os.environ["JSOC_EMAIL"] = "jsoc@sunpy.org"
sphinx_gallery_conf = {
    "backreferences_dir": Path("generated") / "modules",
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": Path("..") / "examples",
    "gallery_dirs": Path("generated") / "gallery",
    "matplotlib_animations": True,
    "default_thumb_file": PNG_ICON,
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("aiapy"),
    "only_warn_on_example_error": True,
}
