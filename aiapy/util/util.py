"""
Miscellaneous utility functions.
"""
import drms
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from sunpy.time import parse_time

from aiapy.util.decorators import validate_channel

__all__ = ['sdo_location', 'telescope_number']


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
    keys = drms.Client().query(
        f'aia.lev1[{(t - 3*u.s).utc.isot}/6s]',
        key='T_OBS, HAEX_OBS, HAEY_OBS, HAEZ_OBS'
    )
    if keys is None or len(keys) == 0:
        raise ValueError('No DRMS records near this time')
    # Linear interpolation between the nearest records within the returned set
    times = Time(list(keys['T_OBS']), scale='utc')
    x = np.interp(t.mjd, times.mjd, keys['HAEX_OBS'])
    y = np.interp(t.mjd, times.mjd, keys['HAEY_OBS'])
    z = np.interp(t.mjd, times.mjd, keys['HAEZ_OBS'])
    return SkyCoord(x=x, y=y, z=z, unit=u.m, representation_type='cartesian',
                    frame='heliocentricmeanecliptic', obstime=t)


@u.quantity_input
@validate_channel('channel')
def telescope_number(channel: u.angstrom):
    """
    For a given channel wavelength, return the associated telescope number.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of AIA channel.

    Returns
    -------
    telescope_number : `int`
        The telescope number of the filter designated by ``channel``
    """
    return {
        94*u.angstrom: 4,
        131*u.angstrom: 1,
        171*u.angstrom: 3,
        193*u.angstrom: 2,
        211*u.angstrom: 2,
        304*u.angstrom: 4,
        335*u.angstrom: 1,
        1600*u.angstrom: 3,
        1700*u.angstrom: 3,
        4500*u.angstrom: 3,
    }[channel]
