"""
Functions for calibrating AIA images.
"""

import copy
import warnings

import numpy as np
from scipy.ndimage import maximum_filter

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_pixel

from sunpy.map import contains_full_disk
from sunpy.map.header_helper import make_fitswcs_header

from aiapy.calibrate.utils import LATEST_CORRECTION_VERSION, _select_epoch_from_correction_table, get_correction_table
from aiapy.utils import AIApyUserWarning, detector_dimensions
from aiapy.utils.decorators import validate_channel

__all__ = ["correct_degradation", "degradation", "register"]


def _target_grid(smap):
    """
    Return the north-aligned, Sun-centered output header and whether ``smap`` fills it.
    """
    detector_shape = tuple(detector_dimensions().to_value(u.pix))
    full_resolution = smap.data.shape == detector_shape
    scale = 0.6 * u.arcsec / u.pix
    binning = 1 if full_resolution else max(1, round((smap.scale.axis1 / scale).to_value(u.one)))
    shape = tuple(np.rint(np.array(detector_shape) / binning).astype(int))
    full_frame = smap.data.shape == shape
    ref_coord = SkyCoord(0, 0, unit="arcsec", frame=smap.coordinate_frame)
    header = make_fitswcs_header(shape, ref_coord, scale=[binning, binning] * scale)
    if not full_frame and contains_full_disk(smap):
        # Full-disk frames resampled by a non-integer factor keep the binned grid when, north-aligned,
        # they overhang it by at most 1% per side, the plate-scale spread full-resolution frames absorb.
        ny, nx = smap.data.shape
        corners = pixel_to_pixel(
            smap.wcs,
            WCS(header),
            np.array([-0.5, nx - 0.5, nx - 0.5, -0.5]),
            np.array([-0.5, -0.5, ny - 0.5, ny - 0.5]),
        )
        full_frame = all(np.abs(c - (n - 1) / 2).max() <= 0.51 * n for c, n in zip(corners, shape[::-1], strict=True))
    return header, full_frame


def register(smap, *, missing=None, order=3, algorithm="interpolation", **kwargs):
    """
    Register a level 1 AIA image or HMI magnetogram using `reproject`.

    Align solar North with the y axis. Full detector frames retain their shape
    with the Sun at the center. Full-resolution pixels are 0.6 arcsec, so HMI
    images land on the same pixel grid as AIA; resampled images use an integer
    multiple of that scale. Full-disk frames resampled by a non-integer factor
    snap to the nearest binned detector grid, centered on the Sun. Cutouts
    share the full-disk grid and expand to contain their transformed edges;
    the Sun need not lie inside them.

    .. note::

        Output headers use standard ``PCi_j`` WCS metadata and may differ
        from the original FITS header.

    Parameters
    ----------
    smap : `~sunpy.map.sources.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`
        A full-disk AIA image or HMI magnetogram, or a cutout of one.
    missing : `float`, optional
        Fill missing pixels with this value, or the minimum non-NaN input
        value if `None`.
    order : `int` or `str`, optional
        Interpolation order, from 0 to 5, defaulting to cubic (3).
        The interpolation names accepted by `reproject` can also be used.
        Ignored for the adaptive and exact algorithms.
    algorithm : `str`, optional
        The reprojection algorithm to use. See the `reproject` documentation
        for details on different algorithms. Default is ``'interpolation'``.
    **kwargs
        Passed to `~sunpy.map.GenericMap.reproject_to`, with
        ``preserve_date_obs=True`` by default. Output timestamps retain the
        input header's time scale.
        Use ``return_footprint=True`` to return coverage.

    Returns
    -------
    registered : `~sunpy.map.sources.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`
        A registered map with processing level 1.5.
    footprint : `~numpy.ndarray`
        Reprojection coverage, returned only when ``return_footprint=True``.
    """
    if smap.processing_level is None or smap.processing_level > 1:
        warnings.warn(
            "Image registration should only be applied to level 1 data",
            AIApyUserWarning,
            stacklevel=2,
        )
    header, full_frame = _target_grid(smap)
    kwargs.setdefault("auto_extent", None if full_frame else "edges")
    kwargs.setdefault("preserve_date_obs", True)
    missing = smap.min() if missing is None else missing
    source = smap
    nan_mask = None
    if algorithm == "interpolation":
        order = {"nearest-neighbor": 0, "bilinear": 1, "biquadratic": 2, "bicubic": 3}.get(order, order)
        kwargs["order"] = order
        if order in (2, 3, 4, 5) and np.isnan(smap.data).any():
            # Spline filtering needs finite input; restore invalid neighborhoods afterward.
            # An order-k spline samples k+1 pixels per axis, so a NaN reaches 2 pixels
            # for orders 2-3 and 3 pixels for orders 4-5.
            nan_mask = maximum_filter(np.isnan(smap.data), size={2: 5, 3: 5, 4: 7, 5: 7}[order], mode="constant")
            source = smap._new_instance(np.where(np.isnan(smap.data), np.nanmedian(smap.data), smap.data), smap.meta)
    return_footprint = kwargs.pop("return_footprint", False)
    result = source.reproject_to(header, algorithm=algorithm, return_footprint=return_footprint, **kwargs)
    smap_l15, footprint = result if return_footprint else (result, None)
    data = smap_l15.data
    if nan_mask is not None:
        mask_map = smap._new_instance(nan_mask.astype(float), smap.meta)
        mask_l15 = mask_map.reproject_to(smap_l15.wcs, order=1, shape_out=data.shape)
        invalid = mask_l15.data > 0
        data[invalid] = np.nan
        if return_footprint:
            footprint[invalid] = 0
    if algorithm == "interpolation":
        np.clip(data, smap.min(), smap.max(), out=data)
    # Fill after clipping to preserve custom missing values.
    data[np.isnan(data)] = missing
    # Reprojection retains only WCS metadata.
    new_meta = copy.deepcopy(smap.meta)
    # Identity matrices are omitted from WCS headers, so remove the old matrix.
    for key in ("crota1", "crota2", "pc1_1", "pc1_2", "pc2_1", "pc2_2", "cd1_1", "cd1_2", "cd2_1", "cd2_2"):
        new_meta.pop(key, None)
    new_meta.update(smap_l15.meta)
    # Reprojection serializes times in UTC; restore the input header's time scale and precision.
    timesys = smap.meta.get("timesys", "UTC").lower()
    date_obs = Time(smap.date if kwargs["preserve_date_obs"] else smap.reference_date, scale=timesys, precision=9)
    date_avg = Time(smap.reference_date, scale=timesys, precision=9)
    new_meta.pop("date_obs", None)
    new_meta.update(
        {
            "date-obs": date_obs.isot,
            "mjd-obs": date_obs.mjd,
            "date-avg": date_avg.isot,
            "mjd-avg": date_avg.mjd,
            "naxis1": data.shape[1],
            "naxis2": data.shape[0],
            "r_sun": (smap.rsun_obs / smap_l15.scale.axis1).to_value(u.pix),
            "lvl_num": 1.5,
            "bitpix": -8 * data.dtype.itemsize,
        }
    )
    newmap = smap._new_instance(data, new_meta, plot_settings=smap.plot_settings)
    return (newmap, footprint) if return_footprint else newmap


def correct_degradation(smap, *, correction_table=None, calibration_version=LATEST_CORRECTION_VERSION):
    """
    Apply time-dependent degradation correction to an AIA map.

    This function applies a time-dependent correction to an AIA observation by
    dividing the observed intensity by the correction factor calculated by
    `degradation`. Any keyword arguments that can be passed to `degradation`
    can also be passed in here.

    Parameters
    ----------
    smap : `~sunpy.map.sources.AIAMap`
        Map to be corrected.
    correction_table : `~astropy.table.Table`
        Table of correction parameters.
        You can get this table by calling `aiapy.calibrate.utils.get_correction_table`.
    calibration_version : `int`, optional
        The version of the correction table to use.
        Defaults to the latest version defined in this module.

    Returns
    -------
    `~sunpy.map.sources.AIAMap`
        Degradation-corrected map.

    See Also
    --------
    degradation
    """
    d = degradation(
        smap.wavelength,
        smap.date,
        correction_table=correction_table,
        calibration_version=calibration_version,
    )
    return smap / d


@u.quantity_input
@validate_channel("channel")
def degradation(
    channel: u.angstrom,
    obstime,
    *,
    correction_table=None,
    calibration_version=LATEST_CORRECTION_VERSION,
) -> u.dimensionless_unscaled:
    r"""
    Correction to account for time-dependent degradation of the instrument.

    The correction factor to account for the time-varying degradation of
    the telescopes is given by a normalization to the calibration epoch
    closest to ``obstime`` and an interpolation within that epoch to
    ``obstime``,

    .. math::

        \frac{A_{eff}(t_{e})}{A_{eff}(t_0)}(1 + p_1\delta t + p_2\delta t^2 + p_3\delta t^3)

    where :math:`A_{eff}(t_e)` is the effective area calculated at the
    calibration epoch for ``obstime``, :math:`A_{eff}(t_0)` is the effective
    area at the first calibration epoch (i.e. at launch),
    :math:`p_1,p_2,p_3` are the interpolation coefficients for the
    ``obstime`` epoch, and :math:`\delta t` is the difference between the
    start time of the epoch and ``obstime``.
    All calibration terms are taken from the ``aia.response`` series in JSOC
    or read from the table input by the user.

    .. note::

        This function is adapted directly from the
        `aia_bp_corrections.pro <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/response/aia_bp_corrections.pro>`__
        routine in SolarSoft.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of the channel.
    obstime : `~astropy.time.Time`
        Observation time.
    correction_table : `~astropy.table.Table`, optional
        Table of correction parameters.
        Defaults to None, which will use the table returned by `aiapy.calibrate.utils.get_correction_table`.
    calibration_version : `int`, optional
        The version of the correction table to use.
        Defaults to the latest version defined in this module.

    Returns
    -------
    `~astropy.units.Quantity`
        Degradation correction factor.

    See Also
    --------
    aiapy.calibrate.utils.get_correction_table
    aiapy.response.Channel.wavelength_response
    aiapy.response.Channel.eve_correction
    """
    if correction_table is None:
        correction_table = get_correction_table()
    if obstime.shape == ():
        obstime = obstime.reshape((1,))
    ratio = np.zeros(obstime.shape)
    poly = np.zeros(obstime.shape)
    for idx, t in enumerate(obstime):
        table = _select_epoch_from_correction_table(
            channel, t, correction_table, calibration_version=calibration_version
        )
        # Time difference between obstime and start of epoch
        dt = (t - table["T_START"][-1]).to(u.day).value
        # Correction to most recent epoch
        ratio[idx] = table["EFF_AREA"][-1] / table["EFF_AREA"][0]
        # Polynomial correction to interpolate within epoch
        poly[idx] = table["EFFA_P1"][-1] * dt + table["EFFA_P2"][-1] * dt**2 + table["EFFA_P3"][-1] * dt**3 + 1.0
    return u.Quantity(poly * ratio)
