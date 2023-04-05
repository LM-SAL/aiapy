from pathlib import Path
from .version import version as __version__

from itertools import compress

_SSW_MIRROR = "https://sohoftp.nascom.nasa.gov/solarsoft/"


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
        ref = textwrap.dedent("\n".join(lines))
    return ref


__citation__ = __bibtex__ = _get_bibtex()
__all__ = ["__version__", "__citation__", "_SSW_MIRROR"]
