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
            attrs.jsoc.Series('aia.response'),
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


def _select_epoch_from_table(channel, obstime, **kwargs):
    """
    Return correction table with only the first epoch and the epoch in
    which `obstime` falls.

    Parameters
    ----------
    obstime: `~astropy.time.Time`
    channel: `~astropy.units.Quantity`
    """
    correction_table = kwargs.get('correction_table', None)
    if isinstance(correction_table, astropy.table.Table):
        table = correction_table
    else:
        table = get_correction_table(correction_table=correction_table)
    # Select only this channel
    table = table[table['WAVE_STR'] == f'{channel.value:.0f}_THIN']
    # Put import here to avoid circular imports
    from aiapy.response.channel import VERSION_NUMBER
    # Select only most recent version number, JSOC keeps some old entries
    table = table[table['VER_NUM'] == VERSION_NUMBER]
    # Select the epoch for the given observation time
    obstime_in_epoch = np.logical_and(obstime >= table['T_START'],
                                      obstime < table['T_STOP'])
    if not obstime_in_epoch.any():
        raise IndexError(f'No valid calibration epoch for {obstime}')
    # Create new table with only first and obstime epochs
    return QTable(table[[0, np.where(obstime_in_epoch)[0][0]]])
