"""
Helper functions for the test suite.
"""

import numpy as np

import astropy.units as u

import sunpy.map

__all__ = ["make_aia_map"]


def make_aia_map(data, wavelength, unit=u.DN / u.s):
    """
    Create a small synthetic AIA map with a minimal, but valid, FITS header.

    Unlike the fixtures in `aiapy.conftest`, which provide realistic level 1
    images, this lets tests control the values, dimensions, and units of the
    input maps, for example to check arithmetic and input validation.
    """
    shape = np.shape(data)
    meta = {
        "telescop": "SDO/AIA",
        "instrume": f"AIA_{wavelength}",
        "detector": "AIA",
        "wavelnth": wavelength,
        "waveunit": "angstrom",
        "bunit": unit.to_string(),
        "exptime": 2.0,
        "crval1": 0.0,
        "crval2": 0.0,
        "crpix1": (shape[1] + 1) / 2,
        "crpix2": (shape[0] + 1) / 2,
        "cdelt1": 0.6,
        "cdelt2": 0.6,
        "ctype1": "HPLN-TAN",
        "ctype2": "HPLT-TAN",
        "cunit1": "arcsec",
        "cunit2": "arcsec",
        "date-obs": "2010-01-01T00:00:00.00",
        "rsun_ref": 696000000.0,
        "rsun_obs": 960.0,
        "dsun_obs": 1.496e11,
        "lvl_num": 1.5,
    }
    return sunpy.map.Map(np.asarray(data, dtype="float32"), meta)
