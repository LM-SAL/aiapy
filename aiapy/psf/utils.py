"""
Contains utility functions for the AIA PSF calculations.
"""

import astropy.units as u

__all__ = ["filter_mesh_parameters"]


def filter_mesh_parameters(*, use_preflightcore=False):
    """
    Geometric parameters for meshes in AIA filters used to calculate the point
    spread function.

    Parameters
    ----------
    use_preflightcore : `bool`, optional
        If True, use the pre-flight values for the filter mesh parameters

    Returns
    -------
    meshinfo : `dict`
        Dictionary with filter mesh information for each channel. Each channel
        entry then contains another dictionary with the following keys
        describing filter mesh properties of that channel
        (see Table 2 of [1]_):

        * ``angle_arm``: Angles of the four entrance filter arms
        * ``error_angle_arm``: Error in angle of the four entrance filter arms
        * ``spacing_e``: Distance between diffraction spikes from entrance filter
        * ``spacing_fp``: Distance between diffraction spikes from focal plane filter
        * ``mesh_pitch``: Pitch of the mesh
        * ``mesh_width``: Width of the mesh
        * ``width``: Width applied to the Gaussian such that *after*
          convolution we have the proper width
          (:math:`4/3` at :math:`1/e` of max)

    References
    ----------
    .. [1] `Grigis, P., Su, Y., Weber M., et al., 2012,
            AIA PSF Characterization and Deconvolution
            <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/DOC/psfreport.pdf>`_

    See Also
    --------
    `calculate_psf` : Calculate the composite point spread function

    Notes
    -----
    These parameters were calculated from the following images and
    reference background images.

    94:

    * image: 'AIA20101016_191039_0094.fits'
    * reference: 'AIA20101016_190903_0094.fits'

    131:

    * image: 'AIA20101016_191035_0131.fits'
    * reference: 'AIA20101016_190911_0131.fits'

    171:

    * image: 'AIA20101016_191037_0171.fits'
    * reference: 'AIA20101016_190901_0171.fits'

    193:

    * image: 'AIA20101016_191056_0193.fits'
    * reference: 'AIA20101016_190844_0193.fits'

    211:

    * image: 'AIA20101016_191038_0211.fits'
    * reference: 'AIA20101016_190902_0211.fits'

    304:

    * image: 'AIA20101016_191021_0304.fits'
    * reference: 'AIA20101016_190845_0304.fits'

    335:

    * image: 'AIA20101016_191041_0335.fits'
    * reference: 'AIA20101016_190905_0335.fits'
    """
    return {
        94 * u.angstrom: {
            "angle_arm": [49.81, 40.16, -40.28, -49.92] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.02, 0.02] * u.degree,
            "spacing_e": 8.99 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.207 * u.pixel,
            "width": (0.951 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.600109, 0.600109] * u.arcsec,
        },
        131 * u.angstrom: {
            "angle_arm": [50.27, 40.17, -39.70, -49.95] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.02, 0.02] * u.degree,
            "spacing_e": 12.37 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.289 * u.pixel,
            "width": (1.033 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.600698, 0.600698] * u.arcsec,
        },
        171 * u.angstrom: {
            "angle_arm": [49.81, 39.57, -40.13, -50.38] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.02, 0.02] * u.degree,
            "spacing_e": 16.26 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.377 * u.pixel,
            "width": (0.962 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.599489, 0.599489] * u.arcsec,
        },
        193 * u.angstrom: {
            "angle_arm": [49.82, 39.57, -40.12, -50.37] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.03, 0.04] * u.degree,
            "spacing_e": 18.39 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.425 * u.pixel,
            "width": (1.512 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.600758, 0.600758] * u.arcsec,
        },
        211 * u.angstrom: {
            "angle_arm": [49.78, 40.08, -40.34, -49.95] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.02, 0.02] * u.degree,
            "spacing_e": 19.97 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.465 * u.pixel,
            "width": (1.199 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.600758, 0.600758] * u.arcsec,
        },
        304 * u.angstrom: {
            "angle_arm": [49.76, 40.18, -40.14, -49.90] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.02, 0.02] * u.degree,
            "spacing_e": 28.87 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.670 * u.pixel,
            "width": (1.247 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.600165, 0.600165] * u.arcsec,
        },
        335 * u.angstrom: {
            "angle_arm": [50.40, 39.80, -39.64, -50.25] * u.degree,
            "error_angle_arm": [0.02, 0.02, 0.02, 0.02] * u.degree,
            "spacing_e": 31.83 * u.pixel,
            "mesh_pitch": 363.0 * u.um,
            "mesh_width": 34.0 * u.um,
            "spacing_fp": 0.738 * u.pixel,
            "width": (0.962 if use_preflightcore else 4.5) * u.pixel,
            "CDELT": [0.600737, 0.600737] * u.arcsec,
        },
    }
