"""Helpers shared by aiapy's tests."""

import numpy as np
from numpy.random import default_rng

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

import sunpy.map
from sunpy.map.header_helper import make_fitswcs_header

from aiapy.utils import telescope_number

RANDOM_GENERATOR = default_rng()
CHANNELS = [94, 131, 171, 193, 211, 304, 335] * u.angstrom


def make_aia_map(
    data,
    wavelength=171,
    *,
    unit=u.DN / u.s,
    level=1.5,
    scale=0.6,
    rotation=0.0,
    exptime=2.0,
    coordinate=None,
    **overrides,
):
    """
    Build a synthetic AIA map with controlled intensities and geometry.

    Channel metadata follows the AIA instrument conventions. ``coordinate`` defaults
    to disk center seen from Earth on 2020-01-01; its ``obstime`` is the exposure
    midpoint and ``DATE-OBS`` is the exposure start. Use real-file fixtures for
    instrument calibration and file-metadata tests. Any other keyword can be replaced
    via ``overrides``.
    """
    data = np.asarray(data)
    if coordinate is None:
        coordinate = SkyCoord(
            0 * u.arcsec, 0 * u.arcsec, frame="helioprojective", observer="earth", obstime="2020-01-01"
        )
    camera = telescope_number(wavelength * u.angstrom)
    header = make_fitswcs_header(
        data.shape,
        coordinate,
        scale=[scale, scale] * u.arcsec / u.pix,
        rotation_angle=rotation * u.deg,
        telescope="SDO/AIA",
        instrument=f"AIA_{camera}",
        wavelength=wavelength * u.angstrom,
    )
    obstime_t = Time(coordinate.obstime.utc, precision=9)
    rsun_obs = header["rsun_obs"] * u.arcsec
    header.update(
        {
            # Instrument and channel
            "detector": "AIA",
            "camera": camera,
            "wave_str": f"{wavelength:g}" + ("" if wavelength in (1600, 1700, 4500) else "_THIN"),
            "waveunit": "angstrom",
            # Observation and calibration
            "exptime": exptime,
            "pixlunit": "DN",
            "bunit": unit.to_string(),
            "lvl_num": level,
            "bitpix": 8 * data.dtype.itemsize * (1 if data.dtype.kind in "iu" else -1),
            "quality": 0,
            "r_sun": (rsun_obs / (scale * u.arcsec / u.pix)).to_value(u.pix),
            "date-obs": (obstime_t - exptime / 2 * u.s).isot,
            "t_obs": f"{obstime_t.isot}Z",
        }
    )
    header.update(overrides)
    return sunpy.map.Map(data, header)
