"""
Miscellaneous utility functions.
"""

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sunpy.time import parse_time

from aiapy.util.decorators import validate_channel
from aiapy.util.net import _get_data_from_jsoc

__all__ = ["check_quality_flag", "sdo_location", "telescope_number"]

# This comes from Table 2 in Section 7.7.6 of the SDO User Guide
_QUALITY_FLAG_MESSAGES = {
    0: "Flatfield data are not available",
    1: "Orbit data are not available",
    2: "Ancillary science data are not available",
    3: "Master pointing data are not available",
    4: "Limb-fit data are not available (generally not applicable to AIA images)",
    8: "Value of MISSVALS keyword is nonzero",
    9: "Value of MISSVALS keyword is more than 1% of the value of the TOTVALS keyword",
    10: "Value of MISSVALS keyword is more than 5% of the value of the TOTVALS keyword",
    11: "Value of MISSVALS keyword is more than 25% of the value of the TOTVALS keyword",
    12: "Spacecraft is not in science pointing mode (coincides with ACS_MODE keyword set to a value other than SCIENCE)",
    13: "Spacecraft eclipse flag is set (coincides with ACS_ECLP keyword set to YES)",
    14: "Spacecraft sun presence flag is not set (coincides with ACS_SUNP keyword set to NO)",
    15: "Spacecraft safe mode flag is set (coincides with ACS_SAFE keyword set to YES)",
    16: "Dark image flag is set (coincides with IMG_TYPE keyword set to DARK)",
    17: "Image Stabilization System (ISS) loop is open (coincides with AISTATE keyword set to OPEN)",
    18: "Calibration image",
    20: "Focus is out of range",
    21: "Register flag is set",
    30: "Quicklook image",
    31: "Image is not available",
}


def check_quality_flag(quality):
    """
    Interpret the ``QUALITY`` flag in the header of an AIA observation.

    AIA images are occasionally affected by operations associated with
    calibration maneuvers or are missing data e.g. due to eclipses. For
    these operating periods, the ``QUALITY`` keyword in the header,
    a 32-bit integer, records the reason for the operation not being
    nominal. Various flags for different operating modes are encoded
    bitwise such that ``QUALITY`` can indicate multiple non-nominal
    operating modes. For more information, see
    `section 7.7.6 of the SDO user guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/>`__.
    This function decodes the ``QUALITY`` flag and returns a string
    indicating the reason(s) for the flags being set. For a nominal
    operation (``QUALITY==0``), "nominal" is returned.

    Parameters
    ----------
    quality: `int`
        Quality flag encoded as an integer. This is typically found in
        the ``QUALITY`` keyword in the FITS header.

    Returns
    -------
    : `list`
        Messages associated with quality flags that are set.
    """
    flags = [i for i in range(32) if (quality & (1 << i))]
    if flags:
        return [_QUALITY_FLAG_MESSAGES.get(f, "(empty)") for f in flags]
    return ["nominal"]


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
    keys = _get_data_from_jsoc(
        query=f"aia.lev1[{(t - 3 * u.s).utc.isot}/6s]", key="T_OBS, HAEX_OBS, HAEY_OBS, HAEZ_OBS"
    )
    if keys is None or len(keys) == 0:
        msg = f"No JSOC records near this time: {t}"
        raise ValueError(msg)
    # Linear interpolation between the nearest records within the returned set
    times = Time(list(keys["T_OBS"]), scale="utc")
    x = np.interp(t.mjd, times.mjd, keys["HAEX_OBS"])
    y = np.interp(t.mjd, times.mjd, keys["HAEY_OBS"])
    z = np.interp(t.mjd, times.mjd, keys["HAEZ_OBS"])
    return SkyCoord(
        x=x,
        y=y,
        z=z,
        unit=u.m,
        representation_type="cartesian",
        frame="heliocentricmeanecliptic",
        obstime=t,
    )


@u.quantity_input
@validate_channel("channel")
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
        94 * u.angstrom: 4,
        131 * u.angstrom: 1,
        171 * u.angstrom: 3,
        193 * u.angstrom: 2,
        211 * u.angstrom: 2,
        304 * u.angstrom: 4,
        335 * u.angstrom: 1,
        1600 * u.angstrom: 3,
        1700 * u.angstrom: 3,
        4500 * u.angstrom: 3,
    }[channel]
