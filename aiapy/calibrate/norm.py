"""
Normalize AIA observations to 1 AU including exposure time correction,
instrument degradation correction, solar disk size unification;
Resize AIA maps and convert the images to Machine-Learning-friendly 
data type (memory map numpy array) with selected meta data information;
"""
from aiapy.calibrate import correct_degradation
from aiapy.calibrate.util import get_correction_table
from sunpy.map.sources.sdo import AIAMap
from sunpy.map import contains_full_disk
from astropy.time import Time
import skimage.transform
import numpy as np
import copy

__all__ = ['norm', 'save_ml']


def norm(smap, disksize=976.0):
    """
    Normalize the AIA map to 1 AU with unified solar disk size.

    This function normalizes the Level 1.5 AIA map to 1 AU including the
    following corrections: exposure time, instrument degradation, and
    disk size variation.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
        A `~sunpy.map.Map` containing a full-disk Level 1.5 AIA image
    disksize : `float`, optional
        Target size of the solar disk in arcsec. The default is set to 976.
    """
    if not isinstance(smap, AIAMap):
        raise ValueError("Input must be an AIAMap.")
    if not contains_full_disk(smap):
        raise ValueError("Input must be a full disk image.")
    if smap.meta['quality'] != 0:
        raise ValueError("Input map has a quality issue.")
    if smap.meta['lvl_num'] < 1.5:
        raise ValueError("Input map needs to be higher than Level 1.0.")

    # correct for degradation
    correction_table = get_correction_table()
    smap = correct_degradation(smap, correction_table=correction_table)

    # get exposure time
    expTime = smap.meta['EXPTIME']

    # get disk size
    rad = smap.meta['RSUN_OBS']
    scale_factor = disksize / rad

    # fix the translation
    X = smap.data
    t = (X.shape[0] / 2.0) - scale_factor * (X.shape[0] / 2.0)

    # normalize by exposure time
    X = X / expTime

    # make a valid mask
    validMask = 1.0 * (X > 0)
    X[np.where(X <= 0.0)] = 0.0

    # rescale the disk size
    XForm = skimage.transform.SimilarityTransform(scale=scale_factor, translation=(t, t))
    Xr = skimage.transform.warp(X, XForm.inverse, preserve_range=True, mode='edge',
                                output_shape=(X.shape[0], X.shape[0]))
    Xd = skimage.transform.warp(validMask, XForm.inverse, preserve_range=True, mode='edge',
                                output_shape=(X.shape[0], X.shape[0]))

    # correct for interpolating over valid pixels
    Xr = np.divide(Xr, (Xd + 1e-8))

    # update meta data
    new_meta = copy.deepcopy(smap.meta)
    new_meta['RSUN_OBS'] = disksize
    new_meta['EXPTIME'] = 1.0

    return smap._new_instance(Xr, new_meta)


def save_ml(smap, outsize=512, filename=''):
    """
    Resize AIA map and convert the image to Machine-Learning-friendly 
    data type (memory map numpy array) with selected meta data information;

    .. note:: The selected meta data information is saved as part of the array:
              FSN, T_OBS (in JD), EXPTIME, WAVELNTH, AECTYPE, NSATPIX (saturation),
              CDELT1, CRPIX1, CDELT2, CRPIX2, CRLN_OBS, CRLT_OBS. For this reason,
              the minimum outsize cannot be smaller than 12.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
        Level 1 or higher AIA image with full frame.
    outsize : `int`, optional
        The output size of the image in pixel, default is set to 512.
    filename : `str`, optional
        The file name of the saved data. If not specified, will be set according
        to the T_OBS of the AIA map: YYYYMMDD_HHMMSS_WAVELNTH.mm

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap`
        Resized AIA map with updated meta data.

    See Also
    --------
    norm
    """
    if outsize < 12:
        raise ValueError("Outsize > 12 is required.")

    X = smap.data

    # figure out the integer factor to downsample by mean
    divideFactor = np.int(np.round(X.shape[0] / outsize))

    # downsample
    Xr = skimage.transform.downscale_local_mean(X, (divideFactor, divideFactor))

    # update meta data
    new_meta = copy.deepcopy(smap.meta)
    new_meta['NAXIS1'] = outsize
    new_meta['NAXIS2'] = outsize
    new_meta['CDELT1'] = smap.meta['CDELT1'] * divideFactor
    new_meta['CRPIX1'] = smap.meta['CRPIX1'] / divideFactor
    new_meta['CDELT2'] = smap.meta['CDELT2'] * divideFactor
    new_meta['CRPIX2'] = smap.meta['CRPIX2'] / divideFactor
    new_meta['X0_MP'] = smap.meta['X0_MP'] / divideFactor
    new_meta['Y0_MP'] = smap.meta['Y0_MP'] / divideFactor

    # if filename is not specified, set to default
    if filename == '':
        obs = new_meta['T_OBS']
        filename = ''.join(('AIA_', obs[:4], obs[5:7], obs[8:10], '_',
                            obs[11:13], obs[14:16], obs[17:19], '_',
                            str(new_meta['wavelnth']), '.mm'))

    # save memory-map to the AIA data array
    newarray = np.memmap(filename, dtype=np.float32, mode='w+', shape=(outsize + 1, outsize))
    newarray[1:, :] = Xr[:, :]

    # save selected AIA map metadata to the first row
    newarray[0, 0] = new_meta['FSN']
    t = Time(new_meta['T_OBS'], format='isot', scale='utc')
    newarray[0, 1] = t.jd
    newarray[0, 2] = new_meta['EXPTIME']
    newarray[0, 3] = new_meta['WAVELNTH']
    newarray[0, 4] = new_meta['AECTYPE']
    newarray[0, 5] = new_meta['NSATPIX']
    newarray[0, 6] = new_meta['CDELT1']
    newarray[0, 7] = new_meta['CRPIX1']
    newarray[0, 8] = new_meta['CDELT2']
    newarray[0, 9] = new_meta['CRPIX2']
    newarray[0, 10] = new_meta['CRLN_OBS']
    newarray[0, 11] = new_meta['CRLT_OBS']
    newarray[0, 12:] = 0.0

    return smap._new_instance(Xr, new_meta)
