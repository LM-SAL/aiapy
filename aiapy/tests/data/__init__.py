from pathlib import Path

from astropy.utils.data import get_pkg_data_filename

__all__ = ["get_test_filepath"]


def get_test_filepath(filename, **kwargs):
    """
    Return the full path to a test file in the ``tests/data`` directory.

    Parameters
    ----------
    filename : `str`
        The name of the file inside the ``data/test`` directory.

    Returns
    -------
    filepath : `str`
        The full path to the file.
    """
    if isinstance(filename, Path):
        # NOTE: get_pkg_data_filename does not accept Path objects
        filename = filename.as_posix()
    return Path(get_pkg_data_filename(filename, package="aiapy.tests.data", **kwargs))
