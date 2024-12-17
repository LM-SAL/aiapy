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

import drms
from sunpy import log
from sunpy.net import attrs as a
from sunpy.net import jsoc

from aiapy import _SSW_MIRRORS
from aiapy.data._manager import manager
from aiapy.util.decorators import validate_channel

__all__ = ["get_correction_table", "get_error_table", "get_pointing_table"]

# Default version of the degradation calibration curve to use.
# This needs to be incremented as the calibration is updated in JSOC.
CALIBRATION_VERSION = 10
# Error table filename available from SSW
AIA_ERROR_FILE = "sdo/aia/response/aia_V{}_error_table.txt"
# Most recent version number for error tables; increment as new versions become available
ERROR_VERSION = 3
# URLs and SHA-256 hashes for each version of the error tables
URL_HASH_ERROR_TABLE = {
    2: (
        [urljoin(mirror, AIA_ERROR_FILE.format(2)) for mirror in _SSW_MIRRORS],
        "ac97ccc48057809723c27e3ef290c7d78ee35791d9054b2188baecfb5c290d0a",
    ),
    3: (
        [urljoin(mirror, AIA_ERROR_FILE.format(3)) for mirror in _SSW_MIRRORS],
        "66ff034923bb0fd1ad20e8f30c7d909e1a80745063957dd6010f81331acaf894",
    ),
}
URL_HASH_POINTING_TABLE = (
    "https://aia.lmsal.com/public/master_aia_pointing3h.csv",
    "a2c80fa0ea3453c62c91f51df045ae04b771d5cbb51c6495ed56de0da2a5482e",
)
URL_HASH_RESPONSE_TABLE = {
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


def get_correction_table(*, correction_table=None, calibration_version=CALIBRATION_VERSION):
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
    correction_table: `str` or `~astropy.table.QTable`, optional
        Path to correction table file or an existing correction table. If None,
        the table will be queried from JSOC. If that fails, the fixed V10
        response table will be used.
    calibration_version : `int`, optional
        The version of the calibration to use when calculating the degradation.
        By default, this is the most recent version available from JSOC. If you
        are using a specific calibration response file, you may need to specify
        this according to the version in that file.
        Defaults to None which will use the most recent version available.

    Returns
    -------
    `~astropy.table.QTable`

    See Also
    --------
    aiapy.calibrate.degradation
    """
    # Fallback for backwards compatibility
    if calibration_version is None:
        calibration_version = CALIBRATION_VERSION
    if correction_table is not None:
        if isinstance(correction_table, QTable):
            return correction_table
        if isinstance(correction_table, str | pathlib.Path):
            table = QTable(astropy_ascii.read(correction_table))
        else:
            msg = "correction_table must be a file path, an existing table, or None."
            raise ValueError(msg)
    else:
        # NOTE: the [!1=1!] disables the drms PrimeKey logic and enables
        # the query to find records that are ordinarily considered
        # identical because the PrimeKeys for this series are WAVE_STR
        # and T_START. Without the !1=1! the query only returns the
        # latest record for each unique combination of those keywords.
        try:
            table = drms.Client().query("aia.response[][!1=1!]", key="**ALL**")
            table = QTable.from_pandas(table)
        except Exception as e:  # NOQA: BLE001
            log.warning("Unable to retrieve response table from JSOC.")
            log.warning(f"Error: {e}")
            log.warning(f"Falling back to fixed V{calibration_version} response table")
            import aiapy.calibrate.util  # NOQA: PLW0406

            table = QTable(
                astropy_ascii.read(getattr(aiapy.calibrate.util, f"fetch_response_table_v{calibration_version}")())
            )
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
def _select_epoch_from_correction_table(channel: u.angstrom, obstime, table, *, version=None):
    """
    Return correction table with only the first epoch and the epoch in which
    ``obstime`` falls and for only one given calibration version.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
    obstime : `~astropy.time.Time`
    table : `~astropy.table.QTable`
    version : `int`
    """
    version = CALIBRATION_VERSION if version is None else version
    # Select only this channel
    # NOTE: The WAVE_STR prime keys for the aia.response JSOC series for the
    # non-EUV channels do not have a thick/thin designation
    thin = "_THIN" if channel not in (1600, 1700, 4500) * u.angstrom else ""
    wave = channel.to(u.angstrom).value
    table = table[table["WAVE_STR"] == f"{wave:.0f}{thin}"]
    table = table[table["VER_NUM"] == version]
    table.sort("DATE")  # Newest entries will be last
    if len(table) == 0:
        extra_msg = " Max version is 3." if channel == 4500 * u.AA else ""
        raise ValueError(
            f"Correction table does not contain calibration for version {version} for {channel}." + extra_msg,
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


def get_pointing_table(start, end):
    """
    Retrieve 3-hourly master pointing table from the JSOC.

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
    end : `~astropy.time.Time`

    Returns
    -------
    `~astropy.table.QTable`

    See Also
    --------
    aiapy.calibrate.update_pointing
    """
    try:
        q = jsoc.JSOCClient().search(
            a.Time(start, end=end),
            a.jsoc.Series.aia_master_pointing3h,
        )
        table = QTable(q)
    except KeyError as e:
        # If there's no pointing information available between these times,
        # JSOC will raise a cryptic KeyError
        # (see https://github.com/LM-SAL/aiapy/issues/71)
        msg = f"Could not find any pointing information between {start} and {end}"
        raise RuntimeError(msg) from e
    except Exception as e:  # NOQA: BLE001
        log.warning("Unable to retrieve pointing table from JSOC.")
        log.warning(f"Error: {e}")
        log.warning("Falling back to fixed pointing table")
        table = QTable(astropy_ascii.read(fetch_pointing_table()))
    table["T_START"] = Time(table["T_START"], scale="utc")
    table["T_STOP"] = Time(table["T_STOP"], scale="utc")
    for c in table.colnames:
        if "X0" in c or "Y0" in c:
            table[c].unit = "pixel"
        if "IMSCALE" in c:
            table[c].unit = "arcsecond / pixel"
        if "INSTROT" in c:
            table[c].unit = "degree"
    # Remove masking on columns with pointing parameters
    for c in table.colnames:
        if any(n in c for n in ["X0", "Y0", "IMSCALE", "INSTROT"]) and hasattr(table[c], "mask"):
            table[c] = table[c].filled(np.nan)
    return table


def get_error_table(error_table=None):
    if error_table is None:
        error_table = fetch_error_table()
    if isinstance(error_table, str | pathlib.Path):
        table = astropy.io.ascii.read(error_table)
    elif isinstance(error_table, QTable):
        table = error_table
    else:
        msg = f"error_table must be a file path, an existing table, or None, not {type(error_table)}"
        raise TypeError(msg)
    table = QTable(table)
    table["DATE"] = Time(table["DATE"], scale="utc")
    table["T_START"] = Time(table["T_START"], scale="utc")
    # NOTE: The warning from erfa here is due to the fact that dates in
    # this table include at least one date from 2030 and converting this
    # date to UTC is ambiguous as the UTC conversion is not well defined
    # at this date.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ErfaWarning)
        table["T_STOP"] = Time(table["T_STOP"], scale="utc")
    table["WAVELNTH"] = u.Quantity(table["WAVELNTH"], "Angstrom")
    table["DNPERPHT"] = u.Quantity(table["DNPERPHT"], "DN photon-1")
    return table


@manager.require("error_table", *URL_HASH_ERROR_TABLE[ERROR_VERSION])
def fetch_error_table():
    return manager.get("error_table")


@manager.require("pointing_table", *URL_HASH_POINTING_TABLE)
def fetch_pointing_table():
    return manager.get("pointing_table")


@manager.require("response_table_v10", *URL_HASH_RESPONSE_TABLE[10])
def fetch_response_table_v10():
    return manager.get("response_table_v10")


@manager.require("response_table_v9", *URL_HASH_RESPONSE_TABLE[9])
def fetch_response_table_v9():
    return manager.get("response_table_v9")


@manager.require("response_table_v8", *URL_HASH_RESPONSE_TABLE[8])
def fetch_response_table_v8():
    return manager.get("response_table_v8")


@manager.require("response_table_v7", *URL_HASH_RESPONSE_TABLE[7])
def fetch_response_table_v7():
    return manager.get("response_table_v7")


@manager.require("response_table_v6", *URL_HASH_RESPONSE_TABLE[6])
def fetch_response_table_v6():
    return manager.get("response_table_v6")


@manager.require("response_table_v4", *URL_HASH_RESPONSE_TABLE[4])
def fetch_response_table_v4():
    return manager.get("response_table_v4")


@manager.require("response_table_v3", *URL_HASH_RESPONSE_TABLE[3])
def fetch_response_table_v3():
    return manager.get("response_table_v3")


@manager.require("response_table_v2", *URL_HASH_RESPONSE_TABLE[2])
def fetch_response_table_v2():
    return manager.get("response_table_v2")
