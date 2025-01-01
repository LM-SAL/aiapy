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

from aiapy import _SSW_MIRRORS
from aiapy.data._manager import manager
from aiapy.util.decorators import validate_channel
from aiapy.util.net import get_data_from_jsoc

__all__ = [
    "get_correction_table",
    "get_error_table",
    "get_pointing_table",
]

# Error table filename available from SSW
_AIA_ERROR_FILE = "sdo/aia/response/aia_V{}_error_table.txt"
# URLs and SHA-256 hashes for each version of the error tables
_URL_HASH_ERROR_TABLE = {
    2: (
        [urljoin(mirror, _AIA_ERROR_FILE.format(2)) for mirror in _SSW_MIRRORS],
        "ac97ccc48057809723c27e3ef290c7d78ee35791d9054b2188baecfb5c290d0a",
    ),
    3: (
        [urljoin(mirror, _AIA_ERROR_FILE.format(3)) for mirror in _SSW_MIRRORS],
        "66ff034923bb0fd1ad20e8f30c7d909e1a80745063957dd6010f81331acaf894",
    ),
}
_URL_HASH_POINTING_TABLE = (
    "https://aia.lmsal.com/public/master_aia_pointing3h.csv",
    "a2c80fa0ea3453c62c91f51df045ae04b771d5cbb51c6495ed56de0da2a5482e",
)
_URL_HASH_RESPONSE_TABLE = {
    10: (
        [urljoin(mirror, "sdo/aia/response/aia_V10_20201119_190000_response_table.txt") for mirror in _SSW_MIRRORS],
        "0a3f2db39d05c44185f6fdeec928089fb55d1ce1e0a805145050c6356cbc6e98",
    ),
    9: (
        [urljoin(mirror, "sdo/aia/response/aia_V9_20200706_215452_response_table.txt") for mirror in _SSW_MIRRORS],
        "f24b384cba9935ae2e8fd3c0644312720cb6add95c49ba46f1961ae4cf0865f9",
    ),
    8: (
        [urljoin(mirror, "sdo/aia/response/aia_V8_20171210_050627_response_table.txt") for mirror in _SSW_MIRRORS],
        "0e8bc6af5a69f80ca9d4fc2a27854681b76574d59eb81d7201b7f618081f0fdd",
    ),
    7: (
        [urljoin(mirror, "sdo/aia/response/aia_V7_20171129_195626_response_table.txt") for mirror in _SSW_MIRRORS],
        "ac2171d549bd6cc6c37e13e505eef1bf0c89fc49bffd037e4ac64f0b895063ac",
    ),
    6: (
        [urljoin(mirror, "sdo/aia/response/aia_V6_20141027_230030_response_table.txt") for mirror in _SSW_MIRRORS],
        "11c148f447d4538db8fd247f74c26b4ae673355e2536f63eb48f9a267e58c7c6",
    ),
    4: (
        [urljoin(mirror, "sdo/aia/response/aia_V4_20130109_204835_response_table.txt") for mirror in _SSW_MIRRORS],
        "7e73f4effa9a8dc55f7b4993a8d181419ef555bf295c4704703ca84d7a0fc3c1",
    ),
    3: (
        [urljoin(mirror, "sdo/aia/response/aia_V3_20120926_201221_response_table.txt") for mirror in _SSW_MIRRORS],
        "0a5d2c2ed1cda18bb9fbdbd51fbf3374e042d20145150632ac95350fc99de68b",
    ),
    2: (
        [urljoin(mirror, "sdo/aia/response/aia_V2_20111129_000000_response_table.txt") for mirror in _SSW_MIRRORS],
        "d55ccd6cb3cb4bd1c688f8663f942f8a872c918a2504e5e474aa97dff45b62c9",
    ),
}


def _fetch_response_table(version: int):
    # Until the delayed feature from sunpy (v6.1) is out, this function
    # will need to be like this.
    if version not in _URL_HASH_RESPONSE_TABLE:
        msg = f"Invalid response table version: {version}"
        raise ValueError(msg)

    @manager.require("response_table_v10", *_URL_HASH_RESPONSE_TABLE[10])
    def fetch_response_table_v10():
        return manager.get("response_table_v10")

    @manager.require("response_table_v9", *_URL_HASH_RESPONSE_TABLE[9])
    def fetch_response_table_v9():
        return manager.get("response_table_v9")

    @manager.require("response_table_v8", *_URL_HASH_RESPONSE_TABLE[8])
    def fetch_response_table_v8():
        return manager.get("response_table_v8")

    @manager.require("response_table_v7", *_URL_HASH_RESPONSE_TABLE[7])
    def fetch_response_table_v7():
        return manager.get("response_table_v7")

    @manager.require("response_table_v6", *_URL_HASH_RESPONSE_TABLE[6])
    def fetch_response_table_v6():
        return manager.get("response_table_v6")

    @manager.require("response_table_v4", *_URL_HASH_RESPONSE_TABLE[4])
    def fetch_response_table_v4():
        return manager.get("response_table_v4")

    @manager.require("response_table_v3", *_URL_HASH_RESPONSE_TABLE[3])
    def fetch_response_table_v3():
        return manager.get("response_table_v3")

    @manager.require("response_table_v2", *_URL_HASH_RESPONSE_TABLE[2])
    def fetch_response_table_v2():
        return manager.get("response_table_v2")

    return locals()[f"fetch_response_table_v{version}"]()


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
    source: pathlib.Path, str or int
        The source of the correction table. If it is a `pathlib.Path`, it must be a file.
        A string file path will error as an invalid source.
        If a string, it must be "jsoc" which will fetch the most recent version from the JSOC,
        otherwise an integer: 3, 4, 6, 7, 8, 9, 10 to use fixed version files from SSW.

    Returns
    -------
    `~astropy.table.QTable`
        Table of degradation correction factors.

    See Also
    --------
    aiapy.calibrate.degradation
    """
    if isinstance(source, pathlib.Path):
        table = QTable(astropy.io.ascii.read(source))
    elif source in _URL_HASH_RESPONSE_TABLE:
        table = QTable(astropy.io.ascii.read(_fetch_response_table(source)))
    elif isinstance(source, str) and source.lower() == "jsoc":
        # NOTE: the [!1=1!] disables the drms PrimeKey logic and enables
        # the query to find records that are ordinarily considered
        # identical because the PrimeKeys for this series are WAVE_STR
        # and T_START. Without the !1=1! the query only returns the
        # latest record for each unique combination of those keywords.
        table = QTable.from_pandas(get_data_from_jsoc(query="aia.response[][!1=1!]", key="**ALL**"))
    else:
        msg = (
            f"correction_table must be a file path (pathlib.Path), 'jsoc' or one of 3, 4, 6, 7, 8, 9, 10. Not {source}"
        )
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
def _select_epoch_from_correction_table(channel: u.angstrom, obstime, table):
    """
    Return correction table with only the first epoch and the epoch in which
    ``obstime`` falls and for only one given calibration version.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
    obstime : `~astropy.time.Time`
    table : `~astropy.table.QTable`
    """
    # Select only this channel
    # NOTE: The WAVE_STR prime keys for the aia.response JSOC series for the
    # non-EUV channels do not have a thick/thin designation
    thin = "_THIN" if channel not in (1600, 1700, 4500) * u.angstrom else ""
    wave = channel.to(u.angstrom).value
    table = table[table["WAVE_STR"] == f"{wave:.0f}{thin}"]
    table.sort("DATE")  # Newest entries will be last
    if len(table) == 0:
        extra_msg = " Max version is 3." if channel == 4500 * u.AA else ""
        raise ValueError(
            f"Correction table does not contain calibration for {channel}." + extra_msg,
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
def fetch_pointing_table():
    return manager.get("pointing_table")


def get_pointing_table(start, end, *, source):
    """
    Retrieve 3-hourly master pointing table from the given source.

    This function can either fetch the pointing table from JSOC or from aia.lmsal.com.

    This function queries `JSOC <http://jsoc.stanford.edu/>`__ for
    the 3-hourly master pointing table (MPT) in the interval defined by
    ``start`` and ``end``.
    The 3-hourly MPT entries are computed from limb fits of images with
    ``T_OBS`` between ``T_START`` and ``T_STOP``.

    .. note:: A MPT entry covers the interval ``[T_START:T_STOP)``;
              that is, the interval includes ``T_START`` and excludes
              ``T_STOP``.

    .. note:: While it is generally true that ``TSTOP = T_START + 3 hours``,
              there are edge cases where ``T_STOP`` is more than 3 hours
              after ``T_START`` because of a calibration, an eclipse,
              or other reasons, but the fits are still calculated based
              on images from ``T_START`` to ``T_START + 3 hours``.
              Pointing is not stable during these periods, so the question
              of which MPT entry to use is not relevant.

    Parameters
    ----------
    start : `~astropy.time.Time`
        Start time of the interval.
    end : `~astropy.time.Time`
        End time of the interval.
    source : str
        Name of the source from which to retrieve the pointing table.
        Must be one of ``"jsoc"`` or ``"lmsal"``.
        Note that the LMSAL pointing table is not updated frequently.

    Returns
    -------
    `~astropy.table.QTable`

    See Also
    --------
    aiapy.calibrate.update_pointing
    """
    if source.lower() == "jsoc":
        table = get_data_from_jsoc(query=f"aia.master_pointing3h[{start.isot}Z-{end.isot}Z]", key="**ALL**")
    elif source.lower() == "lmsal":
        table = QTable(astropy_ascii.read(fetch_pointing_table()))
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


def _fetch_error_table(version: int):
    # Until the delayed feature from sunpy (v6.1) is out, this function
    # will need to be like this.
    @manager.require("error_table_v2", *_URL_HASH_ERROR_TABLE[2])
    def fetch_error_table_v2():
        return manager.get("error_table_v2")

    @manager.require("error_table_v3", *_URL_HASH_ERROR_TABLE[3])
    def fetch_error_table_v3():
        return manager.get("error_table_v3")

    if version == 2:
        return fetch_error_table_v2()
    if version == 3:
        return fetch_error_table_v3()
    msg = f"Invalid error table version: {version}, must be 2 or 3"
    raise ValueError(msg)


def get_error_table(source) -> QTable:
    """
    Fetches the error table from a SSW mirror or uses a local file if one is
    provided.

    Parameters
    ----------
    source : pathlib.Path, int
        If input is a pathlib.Path, it is assumed to be a file path to a local error table.
        If a path is provided as a string, it will error as an invalid source.
        If the input is an int, it is assumed to be a version of the error table (2, 3).

    Returns
    -------
    QTable
        Error table.

    Raises
    ------
    TypeError
        If ``error_table`` is not a file path.
    """
    if isinstance(source, pathlib.Path):
        error_table = QTable(astropy.io.ascii.read(source))
    elif source in [2, 3]:
        error_table = QTable(astropy.io.ascii.read(_fetch_error_table(source)))
    else:
        msg = f"source must be a file path, or  2 or 3, not {source}"
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
