"""
Functions for calibrating AIA images.
"""

import warnings

import numpy as np

import astropy.units as u

from sunpy.map import contains_full_disk
from sunpy.map.sources.sdo import AIAMap, HMIMap
from sunpy.util.decorators import add_common_docstring

from aiapy.calibrate.transform import _rotation_function_names
from aiapy.calibrate.util import _select_epoch_from_correction_table
from aiapy.util import AIApyUserWarning
from aiapy.util.decorators import validate_channel

__all__ = ["correct_degradation", "degradation", "register"]


@add_common_docstring(rotation_function_names=_rotation_function_names)
def register(smap, *, missing=None, order=3, method="scipy"):
    """
    Processes a full-disk level 1 `~sunpy.map.sources.AIAMap` into a level
    1.5 `~sunpy.map.sources.AIAMap`.

    Rotates, scales and translates the image so that solar North is aligned
    with the y axis, each pixel is 0.6 arcsec across, and the center of the
    Sun is at the center of the image. The actual transformation is done by
    the `~sunpy.map.GenericMap.rotate` method.

    .. warning::

        This function might not return a 4096 by 4096 data array
        due to the nature of rotating and scaling the image.
        If you need a 4096 by 4096 image, you will need to pad the array manually,
        update header: crpix1 and crpix2 by the difference divided by 2 in size along that axis.
        Then create a new map.

        Please open an issue on the `aiapy GitHub page <https://github.com/LM-SAL/aiapy/issues>`__
        if you would like to see this changed.

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
    order : `int`, optional
        Order of the spline interpolation.
    method : {{{rotation_function_names}}}, optional
        Rotation function to use. Defaults to ``'scipy'``.

    Returns
    -------
    `~sunpy.map.sources.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`:
        A level 1.5 copy of `~sunpy.map.sources.AIAMap` or
        `~sunpy.map.sources.sdo.HMIMap`.
    """
    # This implementation is taken directly from the `aiaprep` method in
    # sunpy.instr.aia.aiaprep under the terms of the BSD 2 Clause license.
    # See licenses/SUNPY.rst.
    if not isinstance(smap, AIAMap | HMIMap):
        msg = "Input must be an AIAMap or HMIMap."
        raise TypeError(msg)
    if not contains_full_disk(smap):
        msg = "Input must be a full disk image."
        raise ValueError(msg)
    if smap.processing_level is None or smap.processing_level > 1:
        warnings.warn(
            "Image registration should only be applied to level 1 data",
            AIApyUserWarning,
            stacklevel=3,
        )
    # Target scale is 0.6 arcsec/pixel, but this needs to be adjusted if the
    # map has already been rescaled.
    if (smap.scale[0] / 0.6).round() != 1.0 * u.arcsec / u.pix and smap.data.shape != (
        4096,
        4096,
    ):
        scale = (smap.scale[0] / 0.6).round() * 0.6 * u.arcsec
    else:
        scale = 0.6 * u.arcsec  # pragma: no cover # needs a full res image
    scale_factor = smap.scale[0] / scale
    missing = smap.min() if missing is None else missing
    tempmap = smap.rotate(
        recenter=True,
        scale=scale_factor.value,
        order=order,
        missing=missing,
        method=method,
    )
    # Extract center from padded smap.rotate output
    # crpix1 and crpix2 will be equal (recenter=True), as prep does not work with submaps
    center = np.floor(tempmap.meta["crpix1"])
    range_side = (center + np.array([-1, 1]) * smap.data.shape[0] / 2) * u.pix
    newmap = tempmap.submap(
        u.Quantity([range_side[0], range_side[0]]),
        top_right=u.Quantity([range_side[1], range_side[1]]) - 1 * u.pix,
    )
    newmap.meta["r_sun"] = newmap.meta["rsun_obs"] / newmap.meta["cdelt1"]
    newmap.meta["lvl_num"] = 1.5
    newmap.meta["bitpix"] = -64
    return newmap


def correct_degradation(smap, *, correction_table):
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
        You can get this table by calling `aiapy.calibrate.util.get_correction_table`.

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
    correction_table,
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

    .. note:: This function is adapted directly from the
              `aia_bp_corrections.pro <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/response/aia_bp_corrections.pro>`__
              routine in SolarSoft.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of the channel.
    obstime : `~astropy.time.Time`
        Observation time.
    correction_table : `~astropy.table.Table`
        Table of correction parameters.
        You can get this table by calling `aiapy.calibrate.util.get_correction_table`.

    Returns
    -------
    `~astropy.units.Quantity`
        Degradation correction factor.

    See Also
    --------
    aiapy.calibrate.util.get_correction_table
    aiapy.response.Channel.wavelength_response
    aiapy.response.Channel.eve_correction
    """
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
