# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------
from aiapy import __version__
from astropy.utils.exceptions import AstropyDeprecationWarning
from datetime import datetime
from sunpy import log
import aiapy.data.sample  # NOQA
import os
import warnings
from sunpy.util.exceptions import (
    SunpyDeprecationWarning,
    SunpyPendingDeprecationWarning,
)

log.setLevel("DEBUG")
release = __version__
is_development = ".dev" in __version__
project = "aiapy"
copyright = f"{datetime.now().year}, AIA Instrument Team"
author = "AIA Instrument Team"
os.environ["JSOC_EMAIL"] = "nabil.freij@gmail.com"
os.environ["HIDE_PARFIVE_PROGESS"] = "True"
if not is_development:
    warnings.simplefilter("ignore")
warnings.filterwarnings("error", category=SunpyDeprecationWarning)
warnings.filterwarnings("error", category=SunpyPendingDeprecationWarning)
warnings.filterwarnings("error", category=AstropyDeprecationWarning)

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx_automodapi.automodapi",
    "sphinx_automodapi.smart_resolver",
    "sphinx_changelog",
    "sphinx_panels",
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_gallery.gen_gallery",
]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"
default_role = "obj"

# -- Options for intersphinx extension ---------------------------------------
intersphinx_mapping = {
    "python": (
        "https://docs.python.org/3/",
        (None, "http://data.astropy.org/intersphinx/python3.inv"),
    ),
    "numpy": (
        "https://docs.scipy.org/doc/numpy/",
        (None, "http://data.astropy.org/intersphinx/numpy.inv"),
    ),
    "scipy": (
        "https://docs.scipy.org/doc/scipy/reference/",
        (None, "http://data.astropy.org/intersphinx/scipy.inv"),
    ),
    "matplotlib": (
        "https://matplotlib.org/",
        (None, "http://data.astropy.org/intersphinx/matplotlib.inv"),
    ),
    "astropy": ("http://docs.astropy.org/en/stable/", None),
    "sunpy": ("https://docs.sunpy.org/en/stable/", None),
    "skimage": ("https://scikit-image.org/docs/stable/", None),
    "cupy": ("https://docs.cupy.dev/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
from sunpy_sphinx_theme.conf import *  # NOQA

html_theme_options["logo_url"] = "https://aiapy.readthedocs.io/en/stable/"
graphviz_output_format = "svg"
graphviz_dot_args = [
    "-Nfontsize=10",
    "-Nfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Efontsize=10",
    "-Efontname=Helvetica Neue, Helvetica, Arial, sans-serif",
    "-Gfontsize=10",
    "-Gfontname=Helvetica Neue, Helvetica, Arial, sans-serif",
]
# -- Sphinx-gallery ----------------------------------------------------------
sphinx_gallery_conf = {
    "backreferences_dir": os.path.join("generated", "modules"),
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": os.path.join("..", "examples"),
    "gallery_dirs": os.path.join("generated", "gallery"),
    "matplotlib_animations": True,
    "default_thumb_file": "_static/sdo.png",
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("sunpy"),
    "only_warn_on_example_error": True,
}
