"""
Utilities for computing intensity corrections.
"""

import pathlib
import warnings
from urllib.parse import urljoin

import numpy as np
from erfa.core import ErfaWarning

import astropy.io.ascii
import astropy.units as u
from astropy.io import ascii as astropy_ascii
from astropy.table import QTable
from astropy.time import Time

from sunpy import log
from sunpy.time import TimeRange

from aiapy import _SSW_MIRRORS
from aiapy.data._manager import manager
from aiapy.util.decorators import validate_channel
from aiapy.util.net import _get_data_from_jsoc

__all__ = [
    "get_correction_table",
    "get_error_table",
    "get_pointing_table",
]

# Error table filename available from SSW
_AIA_ERROR_FILE = "sdo/aia/response/aia_V{}_error_table.txt"
# URLs and SHA-256 hashes for each version of the error tables
_URL_HASH_ERROR_TABLE = {
    3: (
        [urljoin(mirror, _AIA_ERROR_FILE.format(3)) for mirror in _SSW_MIRRORS],
        "66ff034923bb0fd1ad20e8f30c7d909e1a80745063957dd6010f81331acaf894",
    )
}
_URL_HASH_POINTING_TABLE = (
    "https://aia.lmsal.com/public/master_aia_pointing3h.csv",
    "a2c80fa0ea3453c62c91f51df045ae04b771d5cbb51c6495ed56de0da2a5482e",
)
_URL_HASH_CORRECTION_TABLE = {
    10: (
        [urljoin(mirror, "sdo/aia/response/aia_V10_20201119_190000_response_table.txt") for mirror in _SSW_MIRRORS],
        "0a3f2db39d05c44185f6fdeec928089fb55d1ce1e0a805145050c6356cbc6e98",
    )
}
_URL_HASH_CORRECTION_TABLE["latest"] = _URL_HASH_CORRECTION_TABLE[10]
_URL_HASH_ERROR_TABLE["latest"] = _URL_HASH_ERROR_TABLE[3]


@manager.require("correction_table_latest", *_URL_HASH_CORRECTION_TABLE["latest"])
def _fetch_correction_table_latest():
    return manager.get("correction_table_latest")


@manager.require("error_table_latest", *_URL_HASH_ERROR_TABLE["latest"])
def _fetch_error_table_latest():
    return manager.get("error_table_latest")


def get_correction_table(source):
    """
    Return table of degradation correction factors.

    This function returns a table of parameters for estimating the
    time-dependent degradation of the instrument. By default, this table
    is queried from ``aia.response`` series in
    `JSOC <http://jsoc.stanford.edu/>`__. The correction table can also be read
    from a file by passing a filepath to ``correction_table``. These files are
    typically included in the SDO tree of an SSW installation in
    ``$SSW/sdo/aia/response/`` with filenames like ``aia_V*_response_table.txt``.

    Parameters
    ----------
    source: pathlib.Path, str
        The source of the correction table.
        It can be a string or `pathlib.Path` for a file .
        Otherwise, it must either be "JSOC" which will fetch the most recent version from the JSOC or "SSW" which will fetch the most recent version from SSW.

    Returns
    -------
    `~astropy.table.QTable`
        Table of degradation correction factors.

    See Also
    --------
    aiapy.calibrate.degradation
    """
    if isinstance(source, str) and source.lower() == "ssw":
        table = QTable(astropy.io.ascii.read(_fetch_correction_table_latest()))
    elif isinstance(source, str) and source.lower() == "jsoc":
        # NOTE: the [!1=1!] disables the drms PrimeKey logic and enables
        # the query to find records that are ordinarily considered
        # identical because the PrimeKeys for this series are WAVE_STR
        # and T_START. Without the !1=1! the query only returns the
        # latest record for each unique combination of those keywords.
        table = QTable.from_pandas(_get_data_from_jsoc(query="aia.response[][!1=1!]", key="**ALL**"))
    elif isinstance(source, pathlib.Path | str):
        table = QTable(astropy.io.ascii.read(source))
    else:
        msg = f"correction_table must be a file path (pathlib.Path), 'JSOC' or 'SSW'. Not {source}"
        raise ValueError(msg)
    selected_cols = [
        "DATE",
        "VER_NUM",
        "WAVE_STR",
        "WAVELNTH",
        "T_START",
        "T_STOP",
        "EFFA_P1",
        "EFFA_P2",
        "EFFA_P3",
        "EFF_AREA",
        "EFF_WVLN",
    ]
    table = table[selected_cols]
    table["T_START"] = Time(table["T_START"], scale="utc")
    # NOTE: The warning from erfa here is due to the fact that dates in
    # this table include at least one date from 2030 and converting this
    # date to UTC is ambiguous as the UTC conversion is not well defined
    # at this date.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ErfaWarning)
        table["T_STOP"] = Time(table["T_STOP"], scale="utc")
    table["WAVELNTH"].unit = "Angstrom"
    table["EFF_WVLN"].unit = "Angstrom"
    table["EFF_AREA"].unit = "cm2"
    return table


@u.quantity_input
@validate_channel("channel")
def _select_epoch_from_correction_table(channel: u.angstrom, obstime, correction_table):
    """
    Return correction table with only the first epoch and the epoch in which
    ``obstime`` falls and for only one given calibration version.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
    obstime : `~astropy.time.Time`
    correction_table : `~astropy.table.QTable`
    """
    # Select only this channel
    # NOTE: The WAVE_STR prime keys for the aia.response JSOC series for the
    # non-EUV channels do not have a thick/thin designation
    thin = "_THIN" if channel not in (1600, 1700, 4500) * u.angstrom else ""
    wave = channel.to(u.angstrom).value
    table = correction_table[correction_table["WAVE_STR"] == f"{wave:.0f}{thin}"]
    table.sort("DATE")  # Newest entries will be last
    if len(table) == 0:
        extra_msg = " Max version is 3." if channel == 4500 * u.AA else ""
        msg = f"Correction table does not contain calibration for {channel}.{extra_msg}"
        raise ValueError(
            msg,
        )
    # Select the epoch for the given observation time
    obstime_in_epoch = np.logical_and(obstime >= table["T_START"], obstime < table["T_STOP"])
    if not obstime_in_epoch.any():
        msg = f"No valid calibration epoch for {obstime}"
        raise ValueError(msg)
    # NOTE: In some cases, there may be multiple entries for a single epoch. We want to
    # use the most up-to-date one.
    i_epoch = np.where(obstime_in_epoch)[0]
    if i_epoch.shape[0] > 1:
        log.debug(
            f"Multiple valid epochs for {obstime}. Using the most recent one",
        )
    # Create new table with only first and obstime epochs
    return QTable(table[[0, i_epoch[-1]]])


@manager.require("pointing_table", *_URL_HASH_POINTING_TABLE)
def _fetch_pointing_table():
    return manager.get("pointing_table")


def _get_time(time_range: Time | TimeRange | list | tuple):
    if isinstance(time_range, (Time, list, tuple)):
        start, end = time_range[0], time_range[-1]
    elif isinstance(time_range, TimeRange):
        start, end = time_range.start, time_range.end
    else:
        msg = (
            "time_range must be a Time (with two elements) or TimeRange object, or a tuple of two astropy Time objects"
        )
        raise TypeError(msg)
    return start, end


def get_pointing_table(source, *, time_range=None):
    """
    Retrieve 3-hourly master pointing table from the given source.

    This function can either fetch the pointing table from JSOC or from aia.lmsal.com.

    This function queries `JSOC <http://jsoc.stanford.edu/>`__ for
    the 3-hourly master pointing table (MPT) in the interval defined by
    ``start`` and ``end``.
    The 3-hourly MPT entries are computed from limb fits of images with
    ``T_OBS`` between ``T_START`` and ``T_STOP``.

    .. note::

        A MPT entry covers the interval ``[T_START:T_STOP)``;
        that is, the interval includes ``T_START`` and excludes
        ``T_STOP``.

    .. note::

        While it is generally true that ``TSTOP = T_START + 3 hours``,
        there are edge cases where ``T_STOP`` is more than 3 hours
        after ``T_START`` because of a calibration, an eclipse,
        or other reasons, but the fits are still calculated based
        on images from ``T_START`` to ``T_START + 3 hours``.
        Pointing is not stable during these periods, so the question
        of which MPT entry to use is not relevant.

    .. note::

        The LMSAL pointing table is a static copy of the JSOC table dated from 11/20/2024.
        It was designed as a stop-gap measure while the JSOC recovered from its recent
        water damage.

    Parameters
    ----------
    source : str
        Name of the source from which to retrieve the pointing table.
        Must be one of ``"jsoc"`` or ``"lmsal"``.
    time_range : `~astropy.time.Time`, `~sunpy.time.TimeRange`, optional
        Time range for which to retrieve the pointing table.
        You can pass in a `~astropy.time.Time` object or a tuple of start and end times.
        Alternatively, you can pass in a `~sunpy.time.TimeRange` object.

    Returns
    -------
    `~astropy.table.QTable`

    See Also
    --------
    aiapy.calibrate.update_pointing
    """
    if source.lower() == "jsoc":
        if time_range is None:
            msg = "time_range must be provided if the source is 'jsoc'"
            raise ValueError(msg)
        start, end = _get_time(time_range)
        table = QTable.from_pandas(
            _get_data_from_jsoc(query=f"aia.master_pointing3h[{start.isot}Z-{end.isot}Z]", key="**ALL**")
        )
    elif source.lower() == "lmsal":
        table = QTable(astropy_ascii.read(_fetch_pointing_table()))
    else:
        msg = f"Invalid source: {source}, must be one of 'jsoc' or 'lmsal'"
        raise ValueError(msg)
    table["T_START"] = Time(table["T_START"], scale="utc")
    table["T_STOP"] = Time(table["T_STOP"], scale="utc")
    for c in table.colnames:
        if "X0" in c or "Y0" in c:
            table[c].unit = "pixel"
        if "IMSCALE" in c:
            table[c].unit = "arcsecond / pixel"
        if "INSTROT" in c:
            table[c].unit = "degree"
    for c in table.colnames:
        # Remove masking on columns with pointing parameters
        if any(n in c for n in ["X0", "Y0", "IMSCALE", "INSTROT"]) and hasattr(table[c], "mask"):
            table[c] = table[c].filled(np.nan)
    return table


def get_error_table(source="SSW") -> QTable:
    """
    Fetches the error table from a SSW mirror or uses a local file if one is
    provided.

    Parameters
    ----------
    source : pathlib.Path, str, optional
        The input is allowed to be "SSW" (the default) which will
        fetch the most recent version from SSW.
        Otherwise, it must be a file path.

    Returns
    -------
    `~astropy.table.QTable`
        Error table.

    Raises
    ------
    TypeError
        If ``error_table`` is not a file path.
    """
    if isinstance(source, str) and source.lower() == "ssw":
        error_table = QTable(astropy.io.ascii.read(_fetch_error_table_latest()))
    elif isinstance(source, pathlib.Path | str):
        error_table = QTable(astropy.io.ascii.read(source))
    else:
        msg = f"source must be a filepath, or 'SSW', not {source}"
        raise TypeError(msg)
    error_table["DATE"] = Time(error_table["DATE"], scale="utc")
    error_table["T_START"] = Time(error_table["T_START"], scale="utc")
    # NOTE: The warning from erfa here is due to the fact that dates in
    # this table include at least one date from 2030 and converting this
    # date to UTC is ambiguous as the UTC conversion is not well defined
    # at this date.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ErfaWarning)
        error_table["T_STOP"] = Time(error_table["T_STOP"], scale="utc")
    error_table["WAVELNTH"] = u.Quantity(error_table["WAVELNTH"], "Angstrom")
    error_table["DNPERPHT"] = u.Quantity(error_table["DNPERPHT"], "DN photon-1")
    return error_table
