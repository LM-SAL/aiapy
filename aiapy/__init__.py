"""
``aiapy``
=========

``aiapy`` is a Python package for analyzing data from the Atmospheric Imaging Assembly (AIA) instrument onboard NASA's Solar Dynamics Observatory (SDO) spacecraft.

For more information, see the `aiapy documentation <https://aiapy.readthedocs.io/en/latest/>`__.

For some examples of using aiapy, see `our gallery <https://aiapy.readthedocs.io/en/latest/generated/gallery/index.html>`__.

- `Homepage <https://aia.lmsal.com>`__
"""

from itertools import compress
from pathlib import Path

from .version import version as __version__

_SSW_MIRRORS = [
    # Github mirror managed by LMSAL
    "https://raw.githubusercontent.com/LM-SAL/backup-files/refs/heads/aia/static/",
    "https://soho.nascom.nasa.gov/solarsoft/",
    "https://hesperia.gsfc.nasa.gov/ssw/",
]


def _get_bibtex():
    import textwrap  # NOQA: PLC0415

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
