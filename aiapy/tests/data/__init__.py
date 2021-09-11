from importlib import resources

__all__ = ['get_test_filepath']


def get_test_filepath(filename, **kwargs):
    """
    Return the full path to a test file in the ``data/test`` directory.

    Parameters
    ----------
    filename : `str`
        The name of the file inside the ``data/test`` directory.

    Return
    ------
    filepath : `str`
        The full path to the file.
    """
    with resources.path("aiapy.tests.data", filename) as path:
        return path
