"""
Functions for calibrating AIA images
"""
import copy
import warnings

import numpy as np

import astropy.units as u
from sunpy.map import contains_full_disk
from sunpy.map.sources.sdo import AIAMap, HMIMap

from aiapy.util import AiapyUserWarning
from aiapy.util.decorators import validate_channel
from .util import _select_epoch_from_correction_table, get_correction_table

__all__ = ['register', 'correct_degradation', 'degradation', 'normalize_exposure']


def register(smap, missing=None, order=3, use_scipy=False):
    """
    Processes a full-disk level 1 `~sunpy.map.sources.sdo.AIAMap` into a level
    1.5 `~sunpy.map.sources.sdo.AIAMap`.

    Rotates, scales and translates the image so that solar North is aligned
    with the y axis, each pixel is 0.6 arcsec across, and the center of the
    Sun is at the center of the image. The actual transformation is done by
    the `~sunpy.map.mapbase.GenericMap.rotate` method.


    .. note:: This routine modifies the header information to the standard
              PCi_j WCS formalism. The FITS header resulting in saving a file
              after this procedure will therefore differ from the original
              file.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`
        A `~sunpy.map.Map` containing a full-disk AIA image or HMI magnetogram
    missing : `float`, optional
        If there are missing values after the interpolation, they will be
        filled in with `missing`. If None, the default value will be the
        minimum value of `smap`
    order : `int`, optional
        Order of the spline interpolation
    use_scipy : `bool`, optional
        If True, use `~scipy.ndimage.interpolation.affine_transform` to do the
        image warping. Otherwise, use `~skimage.transform.warp` (recommended).

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`:
        A level 1.5 copy of `~sunpy.map.sources.sdo.AIAMap` or
        `~sunpy.map.sources.sdo.HMIMap`.
    """
    # This implementation is taken directly from the `aiaprep` method in
    # sunpy.instr.aia.aiaprep under the terms of the BSD 2 Clause license.
    # See licenses/SUNPY.rst.
    if not isinstance(smap, (AIAMap, HMIMap)):
        raise ValueError("Input must be an AIAMap or HMIMap.")
    if not contains_full_disk(smap):
        raise ValueError("Input must be a full disk image.")
    # FIXME: this should not be needed. Additional prep calls should
    # not have any effect. This is a precaution in case additional
    # calls to rotate introduce artifacts.
    if smap.processing_level is None or smap.processing_level > 1:
        warnings.warn(
            'Image registration should only be applied to level 1 data',
            AiapyUserWarning
        )
    # Target scale is 0.6 arcsec/pixel, but this needs to be adjusted if the
    # map has already been rescaled.
    if ((smap.scale[0] / 0.6).round() != 1.0 * u.arcsec / u.pix
            and smap.data.shape != (4096, 4096)):
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
        use_scipy=use_scipy
    )
    # extract center from padded smap.rotate output
    # crpix1 and crpix2 will be equal (recenter=True), as prep does not
    # work with submaps
    center = np.floor(tempmap.meta['crpix1'])
    range_side = (center + np.array([-1, 1]) * smap.data.shape[0] / 2) * u.pix
    newmap = tempmap.submap(
        u.Quantity([range_side[0], range_side[0]]),
        top_right=u.Quantity([range_side[1], range_side[1]]) - 1*u.pix)
    newmap.meta['r_sun'] = newmap.meta['rsun_obs'] / newmap.meta['cdelt1']
    newmap.meta['lvl_num'] = 1.5
    newmap.meta['bitpix'] = -64
    return newmap


def correct_degradation(smap, **kwargs):
    """
    Apply time-dependent degradation correction to an AIA map.

    This function applies a time-dependent correction to an AIA observation by
    dividing the observed intensity by the correction factor calculated by
    `degradation`. Any keyword arguments that can be passed to `degradation`
    can also be passed in here.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap`

    See Also
    --------
    degradation
    """
    d = degradation(smap.wavelength, smap.date, **kwargs)
    return smap._new_instance(smap.data / d, smap.meta)


@u.quantity_input
@validate_channel('channel')
def degradation(channel: u.angstrom, obstime,
                **kwargs) -> u.dimensionless_unscaled:
    r"""
    Correction to account for time-dependent degradation of the instrument.

    The correction factor to account for the time-varying degradation of
    the telescopes is given by a normalization to the calibration epoch
    closest to `obstime` and an interpolation within that epoch to
    `obstime`,

    .. math::

        \\frac{A_{eff}(t_{e})}{A_{eff}(t_0)}(1 + p_1\delta t + p_2\delta t^2 + p_3\delta t^3)

    where :math:`A_{eff}(t_e)` is the effective area calculated at the
    calibration epoch for `obstime`, :math:`A_{eff}(t_0)` is the effective
    area at the first calibration epoch (i.e. at launch),
    :math:`p_1,p_2,p_3` are the interpolation coefficients for the
    `obstime` epoch, and :math:`\delta t` is the difference between the
    start time of the epoch and `obstime`.

    All calibration terms are taken from the `aia.response` series in JSOC
    or read from the table input by the user. This function is adapted
    directly from the
    `aia_bp_corrections.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_corrections.pro>`_
    routine in SolarSoft.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
    obstime : `~astropy.time.Time`
    correction_table : `~astropy.table.Table` or `str`, optional
        Table of correction parameters or path to correction table file.
        If not specified, it will be queried from JSOC. See
        `~aiapy.calibrate.util.get_correction_table` for more information.
        If you are processing many images, it is recommended to
        read the correction table once and pass it with this argument to avoid
        multiple redundant network calls.
    calibration_version : `int`, optional
        The version of the calibration to use when calculating the degradation.
        By default, this is the most recent version available from JSOC. If you
        are using a specific calibration response file, you may need to specify
        this according to the version in that file.

    See Also
    --------
    degradation
    aiapy.calibrate.get_correction_table
    aiapy.response.Channel.wavelength_response
    aiapy.response.Channel.eve_correction
    """
    if obstime.shape == ():
        obstime = obstime.reshape((1,))
    ratio = np.zeros(obstime.shape)
    poly = np.zeros(obstime.shape)
    # Do this outside of the loop to avoid repeated queries
    correction_table = get_correction_table(correction_table=kwargs.get('correction_table'))
    for i, t in enumerate(obstime):
        table = _select_epoch_from_correction_table(
            channel, t, correction_table,
            version=kwargs.get('calibration_version')
        )
        # Time difference between obstime and start of epoch
        dt = (t - table['T_START'][-1]).to(u.day).value
        # Correction to most recent epoch
        ratio[i] = table['EFF_AREA'][-1] / table['EFF_AREA'][0]
        # Polynomial correction to interpolate within epoch
        poly[i] = (
            table['EFFA_P1'][-1]*dt
            + table['EFFA_P2'][-1]*dt**2
            + table['EFFA_P3'][-1]*dt**3
            + 1.
        )
    return u.Quantity(poly * ratio)


def normalize_exposure(smap):
    """
    Apply exposure normalization to an AIA map.

    This function applies exposure normalization to an AIA observation by
    dividing the observed intensity by the exposure value extracted from the
    `smap` header.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
    """
    if not isinstance(smap, AIAMap):
        raise ValueError("Input must be an AIAMap")
    if smap.exposure_time <= 0.0 * u.s:
        warnings.warn(
            "Exposure time is less than or equal to 0.0 seconds.",
            AiapyUserWarning
        )
    newmap = smap._new_instance(smap.data / smap.exposure_time.to(u.s).value,
                                copy.deepcopy(smap.meta))
    newmap.meta['exptime'] = 1.0
    newmap.meta['BUNIT'] = 'ct / s'
    return newmap
