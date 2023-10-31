# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------
from aiapy import __version__
from astropy.utils.exceptions import AstropyDeprecationWarning
from datetime import datetime
import os
import warnings
from sunpy.util.exceptions import (
    SunpyDeprecationWarning,
    SunpyPendingDeprecationWarning,
)
from pathlib import Path
from packaging.version import Version

os.environ["JSOC_EMAIL"] = "jsoc@sunpy.org"
os.environ["HIDE_PARFIVE_PROGESS"] = "True"

release = __version__
aiapy_version = Version(__version__)
is_release = not (aiapy_version.is_prerelease or aiapy_version.is_devrelease)
if is_release:
    warnings.simplefilter("ignore")
warnings.filterwarnings("error", category=SunpyDeprecationWarning)
warnings.filterwarnings("error", category=SunpyPendingDeprecationWarning)
warnings.filterwarnings("error", category=AstropyDeprecationWarning)

# -- General configuration ---------------------------------------------------
project = "aiapy"
copyright = f"{datetime.now().year}, AIA Instrument Team"
author = "AIA Instrument Team"
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
    "sphinxext.opengraph",
    "sphinx_design",
    "sphinx_copybutton",
    "hoverxref.extension",
]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = ".rst"
master_doc = "index"
default_role = "obj"
ogp_image = "https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/raw/main/docs/_static/sdo.png"
ogp_use_first_image = True
ogp_description_length = 160
ogp_custom_meta_tags = [
    '<meta property="og:ignore_canonical" content="true" />',
]
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
# -- Sphinx-gallery ----------------------------------------------------------
sphinx_gallery_conf = {
    "backreferences_dir": Path("generated") / "modules",
    "filename_pattern": "^((?!skip_).)*$",
    "examples_dirs": Path("..") / "examples",
    "gallery_dirs": Path("generated") / "gallery",
    "matplotlib_animations": True,
    "default_thumb_file": "_static/sdo.png",
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("sunpy"),
    "only_warn_on_example_error": True,
}

# -- Options for hoverxref -----------------------------------------------------
# adapted from sphinx-hoverxref conf.py
if os.environ.get("READTHEDOCS"):
    # Building on Read the Docs
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
    #
    # roles with std domain
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "ref": "tooltip",
    "term": "tooltip",
}

# -- Options for sphinx-copybutton ---------------------------------------------
# Python Repl + continuation, Bash, ipython and qtconsole + continuation, jupyter-console + continuation
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
