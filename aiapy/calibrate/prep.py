"""
Functions for calibrating AIA images.
"""

import copy
import warnings

import numpy as np

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord

from sunpy.map import contains_full_disk
from sunpy.map.header_helper import make_fitswcs_header
from sunpy.map.sources.sdo import AIAMap, HMIMap
from sunpy.util.decorators import add_common_docstring

from aiapy.calibrate.transform import _rotation_function_names
from aiapy.calibrate.utils import _select_epoch_from_correction_table, get_correction_table
from aiapy.utils import detector_dimensions
from aiapy.utils.decorators import validate_channel
from aiapy.utils.exceptions import AIApyUserWarning

__all__ = ["correct_degradation", "degradation", "register"]


@add_common_docstring(rotation_function_names=_rotation_function_names)
def register(smap, *, missing=None, algorithm="interpolation", **kwargs):
    """
    Rotates, scales and translates the image so that solar North is aligned
    with the y axis, each pixel is 0.6 arcsec across, and the center of the
    Sun is at the center of the image. The actual transformation is done by
    the `reproject` package.

    .. note::

        This routine modifies the header information to the standard
        ``PCi_j`` WCS formalism. The FITS header resulting in saving a file
        after this procedure will therefore differ from the original
        file.

    Parameters
    ----------
    smap : `~sunpy.map.sources.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`
        A `~sunpy.map.Map` containing a full-disk AIA image or HMI magnetogram
    missing : `float`, optional
        If there are missing values after the interpolation, they will be
        filled in with ``missing``. If `None`, the default value will be the
        minimum value of ``smap``
    algorithm : `str`, optional
        The reprojection algorithm to use. See the `reproject` documentation
        for details on different algorithms. Default is 'interpolation'.
    **kwargs
        Additional keyword arguments are passed to the `sunpy.map.GenericMap.reproject_to`.

    Returns
    -------
    `~sunpy.map.sources.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`:
        A promoted copy of `~sunpy.map.sources.AIAMap` or
        `~sunpy.map.sources.sdo.HMIMap`.
    """
    if not isinstance(smap, (AIAMap, HMIMap)):
        warnings.warn("Input is not an AIAMap or an HMIMap", TypeError, stacklevel=2)
    if not contains_full_disk(smap):
        msg = "Input must be a full disk image."
        raise ValueError(msg) from None
    if smap.processing_level is None or smap.processing_level > 1:
        warnings.warn(
            "Image registration should only be applied to level 1 data",
            AIApyUserWarning,
            stacklevel=2,
        )
    # The output WCS is defined in terms of the full-frame image such that
    # the center of the Sun lies at the center of the pixel array.
    shape_full_disk = detector_dimensions()
    ref_pixel_full_disk = (shape_full_disk - 1 * u.pix) / 2
    ref_coord = SkyCoord(0, 0, unit="arcsec", frame=smap.coordinate_frame)
    scale = [0.6, 0.6] * u.arcsec / u.pixel
    # First, construct the full-frame WCS
    header_l15_full_disk = make_fitswcs_header(
        tuple(shape_full_disk.value),
        ref_coord,
        reference_pixel=ref_pixel_full_disk,
        scale=scale,
        rotation_matrix=np.eye(2),
    )
    wcs_l15_full_disk = astropy.wcs.WCS(header_l15_full_disk, preserve_units=True)
    # Use the full-frame WCS directly so disk center is at image center
    kwargs["return_footprint"] = kwargs.get("return_footprint", False)
    # This was selected as the fastest method in local testing
    kwargs["parallel"] = kwargs.get("parallel", True)
    kwargs["block_size"] = kwargs.get("block_size", (1024, 1024))
    smap_l15 = smap.reproject_to(wcs_l15_full_disk, algorithm=algorithm, **kwargs)
    # Fill in missing values
    data = smap_l15.data
    missing = smap.data.min() if missing is None else missing
    data[np.where(np.isnan(data))] = missing
    # Restore metadata (reproject_to only carries over the WCS keywords)
    new_meta = copy.deepcopy(smap.meta)
    # CROTA2 conflicts with the new diagonalized PCi_j
    new_meta.pop("crota2", None)
    # Update the WCS keywords
    new_meta.update(smap_l15.meta)
    # TODO: check if any other keys need to be manually modified?
    new_meta["lvl_num"] = 1.5
    return smap._new_instance(data, new_meta, plot_settings=smap.plot_settings)


def correct_degradation(smap, *, correction_table=None):
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
    )
    return smap / d


@u.quantity_input
@validate_channel("channel")
def degradation(
    channel: u.angstrom,
    obstime,
    *,
    correction_table=None,
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
        table = _select_epoch_from_correction_table(channel, t, correction_table)
        # Time difference between obstime and start of epoch
        dt = (t - table["T_START"][-1]).to(u.day).value
        # Correction to most recent epoch
        ratio[idx] = table["EFF_AREA"][-1] / table["EFF_AREA"][0]
        # Polynomial correction to interpolate within epoch
        poly[idx] = table["EFFA_P1"][-1] * dt + table["EFFA_P2"][-1] * dt**2 + table["EFFA_P3"][-1] * dt**3 + 1.0
    return u.Quantity(poly * ratio)
