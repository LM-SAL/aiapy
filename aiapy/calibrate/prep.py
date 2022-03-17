"""
Functions for calibrating AIA images
"""
import copy
import warnings

import numpy as np

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import pixel_to_pixel
from sunpy.map.header_helper import make_fitswcs_header
from sunpy.map.sources.sdo import AIAMap, HMIMap

from aiapy.util import AiapyUserWarning, detector_dimensions
from aiapy.util.decorators import validate_channel
from .util import _select_epoch_from_correction_table, get_correction_table

__all__ = ['register', 'correct_degradation', 'degradation', 'normalize_exposure']


def register(smap, missing=None, algorithm='interpolation', **kwargs):
    """
    Processes a full-disk level 1 `~sunpy.map.sources.sdo.AIAMap` into a level
    1.5 `~sunpy.map.sources.sdo.AIAMap`.

    Rotates, scales and translates the image so that solar North is aligned
    with the y axis, each pixel is 0.6 arcsec across, and the center of the
    Sun is at the center of the image. The actual transformation is done by
    the `reproject` package.

    .. note:: This routine modifies the header information to the standard
              ``PCi_j`` WCS formalism. The FITS header resulting in saving a file
              after this procedure will therefore differ from the original
              file.

    Additional keyword arguments are passed through to the reprojection function.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`
        A `~sunpy.map.Map` containing a full-disk AIA image or HMI magnetogram
    missing : `float`, optional
        If there are missing values after the interpolation, they will be
        filled in with `missing`. If None, the default value will be the
        minimum value of `smap`
    algorithm : `str`, optional

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap` or `~sunpy.map.sources.sdo.HMIMap`:
        A level 1.5 copy of `~sunpy.map.sources.sdo.AIAMap` or
        `~sunpy.map.sources.sdo.HMIMap`.
    """
    # This check is here because we assume detector dimensions of (4096,4096)
    # and a CDELT of 0.6 arcsec / pix
    if not isinstance(smap, (AIAMap, HMIMap)):
        warnings.warn('Input is not an AIAMap or an HMIMap',
                      AiapyUserWarning)
    # The output WCS is defined in terms of the full-frame image such that
    # the center of the Sun lies at the center of the pixel array.
    shape_full_disk = detector_dimensions()
    ref_pixel_full_disk = (shape_full_disk - 1*u.pix) / 2
    ref_coord = SkyCoord(0, 0, unit='arcsec', frame=smap.coordinate_frame)
    scale = [0.6, 0.6] * u.arcsec / u.pixel
    # First, construct the full-frame WCS
    header_l15_full_disk = make_fitswcs_header(
        tuple(shape_full_disk.value),
        ref_coord,
        reference_pixel=ref_pixel_full_disk,
        scale=scale,
        rotation_matrix=np.eye(2),
    )
    wcs_l15_full_disk = astropy.wcs.WCS(header_l15_full_disk)

    # Find the bottom left corner of the map in the full-frame WCS
    blc_full_disk = pixel_to_pixel(smap.wcs, wcs_l15_full_disk, 0, 0) * u.pix
    # Calculate distance between full-frame center and bottom left corner
    # This is the location of disk center in the aligned WCS of this map
    ref_pixel = ref_pixel_full_disk - blc_full_disk

    # Construct the L1.5 WCS for this map and reproject
    wcs_l15 = astropy.wcs.WCS(make_fitswcs_header(
        smap.data.shape,
        ref_coord,
        reference_pixel=ref_pixel,
        scale=scale,
        rotation_matrix=np.eye(2),
    ))
    smap_l15 = smap.reproject_to(wcs_l15,
                                 return_footprint=False,
                                 algorithm=algorithm,
                                 **kwargs)
    # Fill in missing values
    data = smap_l15.data
    missing = smap.data.min() if missing is None else missing
    data[np.where(np.isnan(data))] = missing

    # Restore metadata (reproject_to only carries over the WCS keywords)
    new_meta = copy.deepcopy(smap.meta)
    # CROTA2 conflicts with the new diagonalized PCi_j
    new_meta.pop('crota2', None)
    # Update the WCS keywords
    new_meta.update(smap_l15.meta)
    # TODO: check if any other keys need to be manually modified?
    new_meta['lvl_num'] = 1.5

    return smap_l15._new_instance(data,
                                  new_meta,
                                  plot_settings=smap_l15.plot_settings)


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
def degradation(channel: u.angstrom, obstime, **kwargs) -> u.dimensionless_unscaled:
    r"""
    Correction to account for time-dependent degradation of the instrument.

    The correction factor to account for the time-varying degradation of
    the telescopes is given by a normalization to the calibration epoch
    closest to `obstime` and an interpolation within that epoch to
    `obstime`,

    .. math::

        \frac{A_{eff}(t_{e})}{A_{eff}(t_0)}(1 + p_1\delta t + p_2\delta t^2 + p_3\delta t^3)

    where :math:`A_{eff}(t_e)` is the effective area calculated at the
    calibration epoch for `obstime`, :math:`A_{eff}(t_0)` is the effective
    area at the first calibration epoch (i.e. at launch),
    :math:`p_1,p_2,p_3` are the interpolation coefficients for the
    `obstime` epoch, and :math:`\delta t` is the difference between the
    start time of the epoch and `obstime`.
    All calibration terms are taken from the `aia.response` series in JSOC
    or read from the table input by the user.

    .. note:: This function is adapted directly from the
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
