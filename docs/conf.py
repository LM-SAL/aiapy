"""
Configuration file for the Sphinx documentation builder.
"""

import os
import sys
import datetime
import warnings
from pathlib import Path

from packaging.version import Version

# -- Read the Docs Specific Configuration --------------------------------------
# This needs to be done before aiapy or sunpy is imported
os.environ["PARFIVE_HIDE_PROGRESS"] = "True"

# -- Check for dependencies ----------------------------------------------------
from sunpy.util import missing_dependencies_by_extra  # NOQA: E402

missing_requirements = missing_dependencies_by_extra("aiapy")["docs"]
if missing_requirements:
    print(  # NOQA: T201
        f"The {' '.join(missing_requirements.keys())} package(s) could not be found and "
        "is needed to build the documentation, please install the 'docs' requirements.",
    )
    sys.exit(1)

# -- Project information -------------------------------------------------------
project = "aiapy"
author = "AIA Instrument Team"
copyright = f"{datetime.datetime.now().year}, {author}"

# Register remote data option with doctest
import doctest  # NOQA: E402

REMOTE_DATA = doctest.register_optionflag("REMOTE_DATA")

from astropy.utils.exceptions import AstropyDeprecationWarning  # NOQA: E402
from matplotlib import MatplotlibDeprecationWarning  # NOQA: E402
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyPendingDeprecationWarning  # NOQA: E402
from sunpy_sphinx_theme import PNG_ICON  # NOQA: E402

from aiapy import __version__  # NOQA: E402

# The full version, including alpha/beta/rc tags
release = __version__
aiapy_version = Version(__version__)
is_release = not (aiapy_version.is_prerelease or aiapy_version.is_devrelease)

warnings.filterwarnings("error", category=SunpyDeprecationWarning)
warnings.filterwarnings("error", category=SunpyPendingDeprecationWarning)
warnings.filterwarnings("error", category=MatplotlibDeprecationWarning)
warnings.filterwarnings("error", category=AstropyDeprecationWarning)

# For the linkcheck
linkcheck_ignore = [
    r"https://doi.org/\d+",
    r"https://element.io/\d+",
    r"https://github.com/\d+",
    r"https://docs.sunpy.org/\d+",
]
linkcheck_anchors = False

# -- General configuration -----------------------------------------------------
# sphinxext-opengraph
ogp_image = "https://raw.githubusercontent.com/sunpy/sunpy-logo/master/generated/sunpy_logo_word.png"
ogp_use_first_image = True
ogp_description_length = 160
ogp_custom_meta_tags = [
    '<meta property="og:ignore_canonical" content="true" />',
]
suppress_warnings = ["app.add_directive"]
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

# -- Options for hoverxref -----------------------------------------------------
if os.environ.get("READTHEDOCS"):
    hoverxref_api_host = "https://readthedocs.org"

    if os.environ.get("PROXIED_API_ENDPOINT"):
        # Use the proxied API endpoint
        # A RTD thing to avoid a CSRF block when docs are using a custom domain
        hoverxref_api_host = "/_"

hoverxref_auto_ref = False
hoverxref_domains = ["py"]
hoverxref_mathjax = True
hoverxref_modal_hover_delay = 500
hoverxref_tooltip_maxwidth = 600  # RTD main window is 696px
hoverxref_intersphinx = list(intersphinx_mapping.keys())
hoverxref_role_types = {
    # Roles within the py domain
    "attr": "tooltip",
    "class": "tooltip",
    "const": "tooltip",
    "data": "tooltip",
    "exc": "tooltip",
    "func": "tooltip",
    "meth": "tooltip",
    "mod": "tooltip",
    "obj": "tooltip",
    # Roles within the std domain
    "confval": "tooltip",
    "hoverxref": "tooltip",
    "ref": "tooltip",  # Would be used by hoverxref_auto_ref if we set it to True
    "term": "tooltip",
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
    # Comes from the theme.
    "default_thumb_file": PNG_ICON,
    "abort_on_example_error": False,
    "plot_gallery": "True",
    "remove_config_comments": True,
    "doc_module": ("aiapy"),
    "only_warn_on_example_error": True,
}

# -- Options for sphinx-copybutton ---------------------------------------------
# Python Repl + continuation, Bash, ipython and qtconsole + continuation, jupyter-console + continuation
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
