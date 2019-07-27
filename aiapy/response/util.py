"""
Utilities for computing response functions
"""
import astropy.time
import astropy.units as u
from sunpy.net import Fido, attrs

__all__ = ['get_correction_table']


def get_correction_table():
    """
    Query JSOC for the time-dependent correction factor parameters
    """
    now = astropy.time.Time.now()
    q = Fido.search(
        attrs.Time(start=now-100*u.year, end=now+100*u.year),  # FIXME: more accurate query?
        attrs.jsoc.Series('aia.response'),
        attrs.jsoc.Keys(['VER_NUM', 'WAVE_STR', 'T_START', 'T_STOP', 'EFFA_P1', 'EFFA_P2',
                         'EFFA_P3', 'EFF_AREA', 'EFF_WVLN']),
    )
    table = list(q)[0].table
    table['T_START'] = astropy.time.Time(table['T_START'], scale='utc')
    table['T_STOP'] = astropy.time.Time(table['T_STOP'], scale='utc')
    return table
