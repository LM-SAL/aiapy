"""
Functions for calibrating AIA images
"""
import numpy as np
import astropy.units as u
from sunpy.map.sources.sdo import AIAMap, HMIMap
from sunpy.map import contains_full_disk

__all__ = ['register']


def register(smap, missing=None, order=3, use_scipy=False):
    """
    Processes a full-disk level 1 `~sunpy.map.sources.sdo.AIAMap` into a level
    1.5 `~sunpy.map.sources.sdo.AIAMap`.

    Rotates, scales and translates the image so that solar North is aligned
    with the y axis, each pixel is 0.6 arcsec across, and the center of the
    Sun is at the center of the image. The actual transformation is done by
    `~sunpy.map.mapbase.GenericMap.rotate` method.


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

    # Target scale is 0.6 arcsec/pixel, but this needs to be adjusted if the
    # map has already been rescaled.
    if ((smap.scale[0] / 0.6).round() != 1.0 * u.arcsec / u.pix
            and smap.data.shape != (4096, 4096)):
        scale = (smap.scale[0] / 0.6).round() * 0.6 * u.arcsec
    else:
        scale = 0.6 * u.arcsec  # pragma: no cover # can't test this because it needs a full res image
    scale_factor = smap.scale[0] / scale

    missing = smap.min() if missing is None else missing

    tempmap = smap.rotate(recenter=True,
                          scale=scale_factor.value,
                          order=order,
                          missing=missing,
                          use_scipy=use_scipy)

    # extract center from padded smap.rotate output
    # crpix1 and crpix2 will be equal (recenter=True), as prep does not
    # work with submaps
    center = np.floor(tempmap.meta['crpix1'])
    range_side = (center + np.array([-1, 1]) * smap.data.shape[0] / 2) * u.pix
    newmap = tempmap.submap(u.Quantity([range_side[0], range_side[0]]),
                            u.Quantity([range_side[1], range_side[1]]))

    newmap.meta['r_sun'] = newmap.meta['rsun_obs'] / newmap.meta['cdelt1']
    newmap.meta['lvl_num'] = 1.5
    newmap.meta['bitpix'] = -64

    return newmap
