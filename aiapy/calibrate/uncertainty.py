"""
Estimate uncertainty on intensities
"""
import numpy as np

import astropy.units as u

from aiapy.util import telescope_number
from aiapy.util.decorators import validate_channel
from .util import get_error_table

__all__ = ['estimate_error']


@u.quantity_input
@validate_channel('channel')
def estimate_error(counts: u.ct / u.pix, channel: u.angstrom, n_sample=1, include_preflight=False,
                   include_eve=False, include_chianti=False, error_table=None, **kwargs) -> u.ct / u.pix:
    """
    Given an observed number of counts estimate the associated errors.

    Calculate the error for a given set of counts for a given channel. This
    calculation includes errors from the shot noise, read noise, quantization, and
    onboard compression. The calculation can also optionally include
    contributions from the photometric calibration and errors in the atomic data.

    .. note:: This function is adapted directly from the
              `aia_bp_corrections.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_estimate_error.pro>`_
              routine in SolarSoft.

    Parameters
    ----------
    counts : `~astropy.units.Quantity`
        Observed counts. These should NOT be divided by the exposure time.
    channel : `~astropy.units.Quantity`
        The corresponding channel for the observed counts.
    n_sample : `int`, optional
        How many measurements (adjacent pixels, or consecutive images) were
        averaged to produce the measured counts.
    include_preflight : `bool`, optional
        Use the preflight photometric calibration. If True, ``include_eve`` must be False.
    include_eve : `bool`, optional
        Use the EVE photometric calibration. If True, ``include_preflight`` must be False.
    include_chianti : `bool`, optional
        If True, include the atomic data errors from CHIANTI in the uncertainty.
    error_table : `~astropy.table.QTable` or path-like, optional
        Error table to use. Can be an existing table or a path to a file. If an error table
        is not specified, the latest version will be downloaded from SolarSoft. Once you've
        downloaded this once, you won't need to download it again unless the remote version
        changes.

    Returns
    -------
    `~astropy.units.Quantity`

    See Also
    --------
    aiapy.calibrate.util.get_error_table
    """
    counts = np.atleast_1d(counts)
    error_table = get_error_table(error_table=error_table)
    error_table = error_table[error_table['WAVELNTH'] == channel]

    # Shot noise
    # NOTE: pixel and photon are "unitless" so we multiply/divide by these
    # units such that the shot noise has the same units as counts
    pix_per_photon = 1 * u.pixel / u.photon  # use this to get units right
    n_photon = counts / error_table['DNPERPHT'] * pix_per_photon
    shot = np.sqrt(n_photon) * error_table['DNPERPHT'] / np.sqrt(n_sample) / pix_per_photon

    # Dark noise
    # NOTE: The dark error of 0.18 is from an analysis of long-term trends in the residual
    # dark error from 2015.
    dark = 0.18 * u.ct / u.pix

    # Read noise
    if kwargs.get('compare_idl', False):
        # The IDL version hardcodes the read noise as 1.15 DN / pixel so we
        # want to use this when comparing against that code and not at any
        # other time. This is why this option is not documented.
        read_noise = 1.15 * u.ct / u.pix
    else:
        # NOTE: This lookup table comes from Table 6 of Boerner et al. (2012)
        # (http://adsabs.harvard.edu/abs/2012SoPh..275...41B). Each
        # channel is attached to one of the four cameras such that the
        # read noise for a given channel depends on what camera it is on.
        read_noise = {
            1: 1.18 * u.ct / u.pix,
            2: 1.20 * u.ct / u.pix,
            3: 1.15 * u.ct / u.pix,
            4: 1.14 * u.ct / u.pix,
        }[telescope_number(channel)]
    read = read_noise / np.sqrt(n_sample)

    # Quantization
    # NOTE: The 1/sqrt(12) factor is the RMS error due to quantization. Under the assumption that the
    # signal is much larger than the least significant bit (LSB), the quantization is not correlated
    # with the signal and has an approximately uniform distribution. The RMS value is thus the
    # standard deviation of this distribution, approximately 1/sqrt(12) LSB.
    # See https://en.wikipedia.org/wiki/Quantization_(signal_processing)
    quant_rms = (1 / np.sqrt(12)) * u.ct / u.pix
    quant = quant_rms / np.sqrt(n_sample)

    # Onboard compression
    compress = shot / error_table['COMPRESS']
    compress[compress < quant_rms] = quant_rms
    compress[counts < 25 * counts.unit] = 0 * counts.unit
    compress /= np.sqrt(n_sample)

    # Photometric calibration
    if include_eve and include_preflight:
        raise ValueError('Cannot include both EVE and pre-flight correction.')
    calib = 0
    if include_eve:
        calib = error_table['EVEERR']
    elif include_preflight:
        calib = error_table['CALERR']
    calib = calib * counts

    # CHIANTI atomic errors
    chianti = error_table['CHIANTI'] if include_chianti else 0
    chianti = chianti * counts

    error_sum = np.sqrt(shot**2 + dark**2 + read**2 + quant**2 + compress**2 + chianti**2 + calib**2)

    return error_sum
