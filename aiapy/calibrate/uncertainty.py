"""
Estimate uncertainty on intensities
"""
from os import error
import pathlib

import astropy.io
from astropy.table import QTable
from astropy.time import Time
import astropy.units as u
import numpy as np
from numpy.core.fromnumeric import compress
from sunpy.data import manager

from aiapy.util.decorators import validate_channel

__all__ = ['estimate_error']

AIA_ERROR_FILE = 'https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/aia_V{}_error_table.txt'  # What changes with version?
VERSION_NUMBER = 3  # Most recent version number for error tables
# URLs and SHA-256 hashes for each version of the error tables
# The URLs are left as a list so that possible mirrors for these files
# can be specified
URL_HASH = {
    2: ((AIA_ERROR_FILE.format(2)),
        'ac97ccc48057809723c27e3ef290c7d78ee35791d9054b2188baecfb5c290d0a'),
    3: ((AIA_ERROR_FILE.format(3)),
        '66ff034923bb0fd1ad20e8f30c7d909e1a80745063957dd6010f81331acaf894'),
}


@u.quantity_input
@validate_channel('channel')
def estimate_error(counts: u.ct / u.pix, channel: u.angstrom, n_sample=1, include_preflight=False,
                   include_eve=False, include_chianti=False, error_table=None) -> u.ct / u.pix:
    """
    Given an observed number of counts estimate the
    uncertainty in the measurement.

    .. note:: This function is adapted directly from the
              `aia_bp_corrections.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_estimate_error.pro>`_
              routine in SolarSoft.

    Parameters
    ----------
    counts
    channel
    n_sample : `int`
        How many measurements (adjacent pixels, or consecutive images) were
        averaged to produce the measured counts.
    include_preflight : `bool`
    include_eve : `bool`
    include_chianti : `bool`
    error_table

    Returns
    -------
    """
    error_table = get_error_table(error_table)
    error_table = error_table[error_table['WAVELNTH'] == channel]

    counts = np.atleast_1d(counts)
    # Shot noise
    # TODO: Figure out what is going on with units here
    pix_per_photon = 1 * u.pixel / u.photon  # use this to get units right
    n_photon = counts / error_table['DNPERPHT'] * pix_per_photon
    shot = np.sqrt(n_photon) * error_table['DNPERPHT'] / np.sqrt(n_sample) / pix_per_photon
    # Dark noise
    # TODO: explain this number?
    dark = 0.18 * u.ct / u.pix
    # Read noise
    # TODO: explain this number?
    read = (1.15 * u.ct / u.pix) * np.sqrt(n_sample) / n_sample
    # Quantization
    # TODO: explain this number?
    quant = (0.288819 * u.ct / u.pix) * np.sqrt(n_sample) / n_sample
    # Onboard compression
    compress = shot / error_table['COMPRESS']
    compress[compress < 0.288819 * compress.unit] = 0.288819 * compress.unit  # TODO: explain this number?
    compress[counts < 25 * counts.unit] = 0 * counts.unit
    compress *= np.sqrt(n_sample) / n_sample
    # Photometric calibration
    if include_eve and include_preflight:
        raise ValueError('Cannot include both EVE and pre-flight correction.')
    if include_eve:
        calib = error_table['EVEERR'] * counts
    elif include_preflight:
        calib = error_table['CALERR'] * counts
    else:
        calib = 0 * counts.unit
    # CHIANTI
    if include_chianti:
        chianti = error_table['CHIANTI'] * counts
    else:
        chianti = 0 * counts.unit

    error_sum = np.sqrt(shot**2 + dark**2 + read**2 + quant**2 + compress**2 + chianti**2 + calib**2)

    return error_sum


def get_error_table(error_table):
    if error_table is None:
        error_table = fetch_error_table()
    if isinstance(error_table, (str, pathlib.Path)):
        table = astropy.io.ascii.read(error_table)
    elif isinstance(error_table, QTable):
        table = error_table
    else:
        raise ValueError('error_table must be a file path, an existing table, or None.')
    table = QTable(table)
    for col in ['DATE', 'T_START', 'T_STOP']:
        table[col] = Time(table[col], scale='utc')
    table['WAVELNTH'] = u.Quantity(table['WAVELNTH'], 'Angstrom')
    table['DNPERPHT'] = u.Quantity(table['DNPERPHT'], 'ct photon-1')
    return table


@manager.require('error_table', *URL_HASH[VERSION_NUMBER])
def fetch_error_table():
    return manager.get('error_table')
