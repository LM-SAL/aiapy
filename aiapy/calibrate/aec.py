"""
Correct the integrated signal of automatic exposure control (AEC) images.
"""

import numpy as np

import astropy.units as u
from astropy.time import Time

__all__ = ["correct_aec"]


@u.quantity_input
def correct_aec(
    table,
    *,
    dn_bias: u.DN = 0 * u.DN,
    exptime_bias: u.s = 0 * u.s,
    saturation_scale=None,
    saturation_shift: u.s = -12 * u.s,
) -> u.DN / u.s:
    """
    Correct the mean count rate of a sequence of AEC-on and AEC-off images.

    During flares, AIA interleaves short AEC-on exposures with long AEC-off
    exposures. Light curves built from the two disagree because of:

    #. An error in the dark pedestal, which matters more for short exposures.
    #. An error in ``EXPTIME`` of a few ms for short exposures, from how the shutter works.
    #. Residual charge from saturated pixels in the previous image, which biases
       the next image high.

    This function corrects ``DATAMEAN / EXPTIME`` for all three.

    .. note::

        This function is adapted from the ``aia_aec_correct.pro`` IDL routine
        by the AIA team.

    .. warning::

        The corrections use empirical coefficients which have to be tuned for
        each observation. The defaults apply no correction except
        ``DATAMEAN / EXPTIME``.

    Parameters
    ----------
    table : `~astropy.table.Table`
        Time-ordered header keywords of one channel, with the columns ``T_OBS``,
        ``AECTYPE``, ``DATAMEAN``, ``EXPTIME`` and ``NSATPIX``.
        Images with ``AECTYPE > 1`` are taken as AEC-on.
        ``DATAMEAN`` can be replaced with any mean in DN, such as the mean of a region.
    dn_bias : `~astropy.units.Quantity`, optional
        Error in the dark pedestal, subtracted from ``DATAMEAN``.
        It can be measured from a dark corner of the images.
    exptime_bias : `~astropy.units.Quantity`, optional
        Error in ``EXPTIME`` of AEC-on images, subtracted from ``EXPTIME``.
    saturation_scale : `float`, optional
        Number of saturated pixels (``NSATPIX``) per 1 DN/s of residual signal.
        If `None`, no saturation correction is applied.
        Frames outside the time range of the other AEC state are extrapolated linearly, as in the IDL routine.
    saturation_shift : `~astropy.units.Quantity`, optional
        Time lag between a saturated image and its residual signal.

    Returns
    -------
    `~astropy.units.Quantity`
        Corrected mean count rate of each image.
    """
    time = Time(table["T_OBS"])
    seconds = (time - time[0]).to_value(u.s)
    aec_on = np.asarray(table["AECTYPE"]) > 1
    exptime = u.Quantity(table["EXPTIME"], u.s)
    exptime = np.where(aec_on, exptime - exptime_bias, exptime)
    rate = (u.Quantity(table["DATAMEAN"], u.DN) - dn_bias) / exptime
    if saturation_scale is not None:
        # Each image is corrected with the NSATPIX of the images in the other AEC state.
        residual = np.asarray(table["NSATPIX"]) / saturation_scale * u.DN / u.s
        shifted = seconds + saturation_shift.to_value(u.s)
        rate[aec_on] -= _interpol(seconds[aec_on], shifted[~aec_on], residual[~aec_on])
        rate[~aec_on] -= _interpol(seconds[~aec_on], shifted[aec_on], residual[aec_on])
    return rate


def _interpol(x, xp, fp):
    """
    Linear interpolation that extrapolates from the end segments, like IDL's ``INTERPOL``.
    """
    i = np.clip(np.searchsorted(xp, x, side="right") - 1, 0, len(xp) - 2)
    return fp[i] + (x - xp[i]) * (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])
