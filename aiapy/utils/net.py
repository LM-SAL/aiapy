"""
Provides basic functionality for querying the JSOC for internal use.
"""

import drms

__all__ = ["_get_data_from_jsoc"]


def _get_data_from_jsoc(query, *, key, seg=None):
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
        jsoc_result = drms.Client().query(query, key=key, seg=seg)
    except Exception as e:
        msg = f"Unable to query the JSOC.\n Error message: {e}"
        raise OSError(msg) from e
    if len(jsoc_result) == 0:
        msg = f"No data found for this query: {query}, key: {key}, seg: {seg}"
        raise RuntimeError(msg)
    return jsoc_result
