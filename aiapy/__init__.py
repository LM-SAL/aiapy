"""
``aiapy``
=========

Python library for AIA data analysis.

* Homepage: https://aia.lmsal.com
* Documentation: https://aiapy.readthedocs.io/en/stable/
"""

from itertools import compress
from pathlib import Path

from .version import version as __version__

_SSW_MIRRORS = [
    "https://soho.nascom.nasa.gov/solarsoft/",
    "https://hesperia.gsfc.nasa.gov/ssw/",
]


def _get_bibtex():
    import textwrap

    # Set the bibtex entry to the article referenced in CITATION.rst
    citation_file = Path(__file__).parent / "CITATION.rst"

    # Explicitly specify UTF-8 encoding in case the system's default encoding is problematic
    with Path.open(citation_file, "r", encoding="utf-8") as citation:
        # Extract the first bibtex block:
        ref = citation.read().partition(".. code:: bibtex\n\n")[2]
        lines = ref.split("\n")
        # Only read the lines which are indented
        lines = list(compress(lines, [line.startswith("   ") for line in lines]))
        return textwrap.dedent("\n".join(lines))


__citation__ = __bibtex__ = _get_bibtex()
__all__ = ["_SSW_MIRRORS", "__citation__", "__version__"]
