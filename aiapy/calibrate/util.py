"""
Utilities for computing intensity corrections
"""
import numpy as np
from astropy.time import Time
import astropy.units as u
import astropy.io.ascii
from astropy.table import QTable
from sunpy.net import jsoc, attrs

__all__ = ['get_correction_table']

# Default version of the degradation calibration curve to use. This
# needs to be incremented as the calibration is updated in JSOC.
CALIBRATION_VERSION = 9


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
    correction_table: `str`, optional
        Path to correction table file

    See Also
    --------
    aiapy.calibrate.degradation
    """
    if correction_table is not None:
        table = astropy.io.ascii.read(correction_table)
    else:
        now = Time.now()
        q = jsoc.JSOCClient().search_metadata(
            # FIXME: more accurate time range?
            attrs.Time(start=now-100*u.year, end=now+100*u.year),
            # NOTE: the [!1=1!] disables the drms PrimeKey logic and enables
            # the query to find records that are ordinarily considered
            # identical (because the PrimeKeys for this series are WAVE_STR
            # and T_START, so without the !1=1! the query only returns the
            # latest record for each unique combination of those keywords
            attrs.jsoc.Series('aia.response[!1=1!]'),
            attrs.jsoc.Keys(['VER_NUM', 'WAVE_STR', 'T_START', 'T_STOP',
                             'EFFA_P1', 'EFFA_P2', 'EFFA_P3', 'EFF_AREA',
                             'EFF_WVLN']),
        )
        table = QTable.from_pandas(q)

    table['T_START'] = Time(table['T_START'], scale='utc')
    table['T_STOP'] = Time(table['T_STOP'], scale='utc')
    table['WAVELNTH'].unit = 'Angstrom'
    table['EFF_WVLN'].unit = 'Angstrom'
    table['EFF_AREA'].unit = 'cm2'

    return table


@u.quantity_input
def _select_epoch_from_table(channel: u.angstrom, obstime, **kwargs):
    """
    Return correction table with only the first epoch and the epoch in
    which `obstime` falls and for only one given calibration version.

    Parameters
    ----------
    channel: `~astropy.units.Quantity`
    obstime: `~astropy.time.Time`
    """
    correction_table = kwargs.get('correction_table', None)
    if isinstance(correction_table, astropy.table.Table):
        table = correction_table
    else:
        table = get_correction_table(correction_table=correction_table)
    # Select only this channel
    # NOTE: The WAVE_STR prime keys for the aia.response JSOC series for the
    # non-EUV channels do not have a thick/thin designation
    thin = '_THIN' if channel not in (1600, 1700, 4500)*u.angstrom else ''
    wave = channel.to(u.angstrom).value
    table = table[table['WAVE_STR'] == f'{wave:.0f}{thin}']
    if len(table) == 0:
        raise IndexError(f'Correction table does not contain calibration for wavelength {wave:.0f}')
    version = kwargs.get('calibration_version')
    version = CALIBRATION_VERSION if version is None else version
    table = table[table['VER_NUM'] == version]
    if len(table) == 0:
        raise IndexError(f'Correction table does not contain calibration for version {version}')
    # Select the epoch for the given observation time
    obstime_in_epoch = np.logical_and(obstime >= table['T_START'],
                                      obstime < table['T_STOP'])
    if not obstime_in_epoch.any():
        raise IndexError(f'No valid calibration epoch for {obstime}')
    # Create new table with only first and obstime epochs
    return QTable(table[[0, np.where(obstime_in_epoch)[0][0]]])
