"""
Provides basic functionality for querying the JSOC for internal use.
"""

import drms

__all__ = ["get_data_from_jsoc"]


def get_data_from_jsoc(query, *, key, seg=None):
    """
    A simple wrapper around `~drms.Client.query` that raises a more informative
    error message.

    This was created to avoid having to write a try/except block for every call to the JSOC within the package.

    Please note that this function is not intended to be used outside of the package.

    Parameters
    ----------
    query : str
        Query string to be passed to `~drms.Client.query`.
    key : str
        Key to be passed to `~drms.Client.query`.
    seg : str
        Seg to be passed to `~drms.Client.query`.

    Returns
    -------
    pandas.DataFrame
        The results of the query.

    Raises
    ------
    OSError
        If the query fails for any reason.
    """
    try:
        return drms.Client().query(query, key=key, seg=seg)
    except KeyError as e:
        # This probably should not be here but yolo.
        # If there's no pointing information available between these times,
        # JSOC will raise a cryptic KeyError
        # (see https://github.com/LM-SAL/aiapy/issues/71)
        msg = "Could not find any pointing information"
        raise RuntimeError(msg) from e
    except Exception as e:
        msg = f"Unable to query the JSOC.\n Error message: {e}"
        raise OSError(msg) from e
