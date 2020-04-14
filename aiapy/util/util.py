"""
Miscellaneous utility functions
"""
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import drms
import numpy as np
from sunpy.time import parse_time

__all__ = ['sdo_location']


def sdo_location(time):
    """
    Returns the location of SDO at the given time.

    Parameters
    ----------
    time : `astropy.time.Time`
        The time can also be supplied in a format compatible with :func:`sunpy.time.parse_time`.

    Returns
    -------
    `astropy.coordinates.SkyCoord`

    Notes
    -----
    This function performs a linear interpolation between the nearest DRMS records on either side of
    the given time.
    """
    t = parse_time(time)

    # Query for +/- 3 seconds around the given time
    keys = drms.Client().query(f'aia.lev1[{(t - 3*u.s).utc.isot}/6s]',
                               key='T_OBS, HAEX_OBS, HAEY_OBS, HAEZ_OBS')
    if len(keys) == 0:
        raise ValueError('No DRMS records near this time')

    # Linear interpolation between the nearest records within the returned set
    times = Time(list(keys['T_OBS']), scale='utc')
    x = np.interp(t.mjd, times.mjd, keys['HAEX_OBS'])
    y = np.interp(t.mjd, times.mjd, keys['HAEY_OBS'])
    z = np.interp(t.mjd, times.mjd, keys['HAEZ_OBS'])

    return SkyCoord(x=x, y=y, z=z, unit=u.m, representation_type='cartesian',
                    frame='heliocentricmeanecliptic', obstime=t)
