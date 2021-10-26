"""
Utilities for computing intensity corrections.
"""
import pathlib
import warnings

import numpy as np

import astropy.io.ascii
import astropy.units as u
from astropy.table import QTable
from astropy.time import Time
from sunpy.data import manager
from sunpy.net import attrs, jsoc

from aiapy.util.decorators import validate_channel
from aiapy.util.exceptions import AiapyUserWarning

__all__ = ['get_correction_table', 'get_pointing_table', 'get_error_table']

# Default version of the degradation calibration curve to use.
# This needs to be incremented as the calibration is updated in JSOC.
CALIBRATION_VERSION = 10
# Error table filename available from SSW
AIA_ERROR_FILE = 'https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/aia_V{}_error_table.txt'
# Most recent version number for error tables; increment as new versions become available
ERROR_VERSION = 3
# URLs and SHA-256 hashes for each version of the error tables
# The URLs are left as a list so that possible mirrors for these files
# can be specified
URL_HASH = {
    2: ((AIA_ERROR_FILE.format(2)),
        'ac97ccc48057809723c27e3ef290c7d78ee35791d9054b2188baecfb5c290d0a'),
    3: ((AIA_ERROR_FILE.format(3)),
        '66ff034923bb0fd1ad20e8f30c7d909e1a80745063957dd6010f81331acaf894'),
}


def get_correction_table(correction_table=None):
    """
    Return table of degradation correction factors.

    This function returns a table of parameters for estimating the
    time-dependent degradation of the instrument. By default, this table
    is queried from `aia.response` series in
    `JSOC <http://jsoc.stanford.edu/>`_. The correction table can also be read
    from a file by passing a filepath to `correction_table`. These files are
    typically included in the SDO tree of an SSW installation in
    `$SSW/sdo/aia/response/` with filenames like `aia_V*_response_table.txt`.

    Parameters
    ----------
    correction_table: `str` or `~astropy.table.QTable`, optional
        Path to correction table file or an existing correction table. If None,
        the table will be queried from JSOC.

    Returns
    -------
    `~astropy.table.QTable`

    See Also
    --------
    aiapy.calibrate.degradation
    """
    if isinstance(correction_table, astropy.table.QTable):
        return correction_table
    if correction_table is not None:
        if isinstance(correction_table, (str, pathlib.Path)):
            table = QTable(astropy.io.ascii.read(correction_table))
        else:
            raise ValueError('correction_table must be a file path, an existing table, or None.')
    else:
        now = Time.now()
        q = jsoc.JSOCClient().search(
            # FIXME: more accurate time range?
            attrs.Time(start=now-100*u.year, end=now+100*u.year),
            # NOTE: the [!1=1!] disables the drms PrimeKey logic and enables
            # the query to find records that are ordinarily considered
            # identical because the PrimeKeys for this series are WAVE_STR
            # and T_START. Without the !1=1! the query only returns the
            # latest record for each unique combination of those keywords.
            attrs.jsoc.Series('aia.response[!1=1!]'),
        )
        table = q.show(
            'DATE',
            'VER_NUM',
            'WAVE_STR',
            'WAVELNTH',
            'T_START',
            'T_STOP',
            'EFFA_P1',
            'EFFA_P2',
            'EFFA_P3',
            'EFF_AREA',
            'EFF_WVLN',
        )
    table['T_START'] = Time(table['T_START'], scale='utc')
    table['T_STOP'] = Time(table['T_STOP'], scale='utc')
    table['WAVELNTH'].unit = 'Angstrom'
    table['EFF_WVLN'].unit = 'Angstrom'
    table['EFF_AREA'].unit = 'cm2'
    return table


@u.quantity_input
@validate_channel('channel')
def _select_epoch_from_correction_table(channel: u.angstrom, obstime, table, version=None):
    """
    Return correction table with only the first epoch and the epoch in
    which `obstime` falls and for only one given calibration version.

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
    thin = '_THIN' if channel not in (1600, 1700, 4500)*u.angstrom else ''
    wave = channel.to(u.angstrom).value
    table = table[table['WAVE_STR'] == f'{wave:.0f}{thin}']
    table = table[table['VER_NUM'] == version]
    table.sort('DATE')  # Newest entries will be last
    if len(table) == 0:
        raise IndexError(f'Correction table does not contain calibration for version {version}')
    # Select the epoch for the given observation time
    obstime_in_epoch = np.logical_and(obstime >= table['T_START'], obstime < table['T_STOP'])
    if not obstime_in_epoch.any():
        raise IndexError(f'No valid calibration epoch for {obstime}')
    # NOTE: In some cases, there may be multiple entries for a single epoch. We want to
    # use the most up-to-date one.
    i_epoch = np.where(obstime_in_epoch)[0]
    if i_epoch.shape[0] > 1:
        warnings.warn(f'Multiple valid epochs for {obstime}. Using the most recent one', AiapyUserWarning)
    # Create new table with only first and obstime epochs
    return QTable(table[[0, i_epoch[-1]]])


def get_pointing_table(start, end):
    """
    Query JSOC for the most up-to-date pointing information.

    This function queries `JSOC <http://jsoc.stanford.edu/>`__ for
    the 3-hourly pointing information in the interval defined by
    `start` and `end`.

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
    q = jsoc.JSOCClient().search(
        attrs.Time(start, end=end),
        attrs.jsoc.Series.aia_master_pointing3h,
    )
    table = q.show()
    if len(table.columns) == 0:
        # If there's no pointing information available between these times,
        # JSOC will raise a cryptic KeyError
        # (see https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/issues/84)
        raise RuntimeError(
            f'Could not find any pointing information between {start} and {end}'
        )
    table['T_START'] = Time(table['T_START'], scale='utc')
    table['T_STOP'] = Time(table['T_STOP'], scale='utc')
    for c in table.colnames:
        if 'X0' in c or 'Y0' in c:
            table[c].unit = 'arcsecond'
        if 'IMSCALE' in c:
            table[c].unit = 'arcsecond / pixel'
        if 'INSTROT' in c:
            table[c].unit = 'degree'
    return table


def get_error_table(error_table=None):
    if error_table is None:
        error_table = fetch_error_table()
    if isinstance(error_table, (str, pathlib.Path)):
        table = astropy.io.ascii.read(error_table)
    elif isinstance(error_table, QTable):
        table = error_table
    else:
        raise ValueError('error_table must be a file path, an existing table, or None.')
    table = QTable(table)
    for col in ['DATE', 'T_START', 'T_STOP']:
        table[col] = Time(table[col], scale='utc')
    table['WAVELNTH'] = u.Quantity(table['WAVELNTH'], 'Angstrom')
    table['DNPERPHT'] = u.Quantity(table['DNPERPHT'], 'ct photon-1')
    return table


@manager.require('error_table', *URL_HASH[ERROR_VERSION])
def fetch_error_table():
    return manager.get('error_table')
