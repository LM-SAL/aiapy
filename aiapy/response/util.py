"""
Utilities for computing response functions
"""
from astropy.time import Time
import astropy.units as u
from astropy.table import Table
from sunpy.net import jsoc, attrs

__all__ = ['get_correction_table']


def get_correction_table():
    """
    Query JSOC for the time-dependent correction factor parameters
    """
    now = Time.now()
    q = jsoc.JSOCClient().search_metadata(
        # FIXME: more accurate time range?
        attrs.Time(start=now-100*u.year, end=now+100*u.year),
        attrs.jsoc.Series('aia.response'),
        attrs.jsoc.Keys(['VER_NUM', 'WAVE_STR', 'T_START', 'T_STOP', 'EFFA_P1',
                         'EFFA_P2', 'EFFA_P3', 'EFF_AREA', 'EFF_WVLN']),
    )
    table = Table.from_pandas(q)
    table['T_START'] = Time(table['T_START'], scale='utc')
    table['T_STOP'] = Time(table['T_STOP'], scale='utc')
    return table
