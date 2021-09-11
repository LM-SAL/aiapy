"""
Calculate the point spread function (PSF) for the AIA telescopes.
"""
import numpy as np

import astropy.units as u

from aiapy.util.decorators import validate_channel

try:
    import cupy
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

__all__ = ['psf', 'filter_mesh_parameters', '_psf']


def filter_mesh_parameters(use_preflightcore=False):
    """
    Geometric parameters for meshes in AIA filters used to calculate the
    point spread function.

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

        * `angle_arm`: Angles of the four entrance filter arms
        * `error_angle_arm`: Error in angle of the four entrance filter arms
        * `spacing_e`: Distance between diffraction spikes from entrance filter
        * `spacing_fp`: Distance between diffraction spikes from focal plane filter
        * `mesh_pitch`: Pitch of the mesh
        * `mesh_width`: Width of the mesh
        * `width`: Width applied to the Gaussian such that *after*
          convolution we have the proper width
          (:math:`4/3` at :math:`1/e` of max)

    References
    ----------
    .. [1] `Grigis, P., Su, Y., Weber M., et al., 2012,
            AIA PSF Characterization and Deconvolution
            <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/psf/DOC/psfreport.pdf>`_

    See Also
    --------
    psf : Calculate the composite point spread function
    """
    # These parameters were calculated from the following images and
    # reference background images:
    #   94:
    #       image: 'AIA20101016_191039_0094.fits'
    #       reference: 'AIA20101016_190903_0094.fits'
    #   131:
    #       image: 'AIA20101016_191035_0131.fits'
    #       reference: 'AIA20101016_190911_0131.fits'
    #   171:
    #       image: 'AIA20101016_191037_0171.fits'
    #       reference: 'AIA20101016_190901_0171.fits'
    #   193:
    #       image: 'AIA20101016_191056_0193.fits'
    #       reference: 'AIA20101016_190844_0193.fits'
    #   211:
    #       image: 'AIA20101016_191038_0211.fits'
    #       reference: 'AIA20101016_190902_0211.fits'
    #   304:
    #       image: 'AIA20101016_191021_0304.fits'
    #       reference: 'AIA20101016_190845_0304.fits'
    #   335:
    #       image: 'AIA20101016_191041_0335.fits'
    #       reference: 'AIA20101016_190905_0335.fits'
    # TODO: put this in another file, either JSON or asdf
    return {
        94 * u.angstrom: {
            'angle_arm': [49.81, 40.16, -40.28, -49.92] * u.deg,
            'error_angle_arm': [0.02, 0.02, 0.02, 0.02] * u.deg,
            'spacing_e': 8.99 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.207 * u.pixel,
            'width': (0.951 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.600109, 0.600109]*u.arcsec,
        },
        131 * u.angstrom: {
            'angle_arm': [50.27, 40.17, -39.70, -49.95] * u.deg,
            'error_angle_arm': [0.02, 0.02, 0.02, 0.02] * u.deg,
            'spacing_e': 12.37 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.289 * u.pixel,
            'width': (1.033 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.600698, 0.600698]*u.arcsec,
        },
        171 * u.angstrom: {
            'angle_arm': [49.81, 39.57, -40.13, -50.38] * u.deg,
            'error_angle_arm': [0.02, 0.02, 0.02, 0.02] * u.deg,
            'spacing_e': 16.26 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.377 * u.pixel,
            'width': (0.962 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.599489, 0.599489]*u.arcsec,
        },
        193 * u.angstrom: {
            'angle_arm': [49.82, 39.57, -40.12, -50.37] * u.deg,
            'error_angle_arm': [0.02, 0.02, 0.03, 0.04] * u.deg,
            'spacing_e': 18.39 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.425 * u.pixel,
            'width': (1.512 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.600758, 0.600758]*u.arcsec,
        },
        211 * u.angstrom: {
            'angle_arm': [49.78, 40.08, -40.34, -49.95] * u.deg,
            'error_angle_arm': [0.02, 0.02, 0.02, 0.02] * u.deg,
            'spacing_e': 19.97 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.465 * u.pixel,
            'width': (1.199 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.600758, 0.600758]*u.arcsec,
        },
        304 * u.angstrom: {
            'angle_arm': [49.76, 40.18, -40.14, -49.90] * u.degree,
            'error_angle_arm': [0.02, 0.02, 0.02, 0.02] * u.deg,
            'spacing_e': 28.87 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.670 * u.pixel,
            'width': (1.247 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.600165, 0.600165]*u.arcsec,
        },
        335 * u.angstrom: {
            'angle_arm': [50.40, 39.80, -39.64, -50.25] * u.degree,
            'error_angle_arm': [0.02, 0.02, 0.02, 0.02] * u.deg,
            'spacing_e': 31.83 * u.pixel,
            'mesh_pitch': 363.0 * u.um,
            'mesh_width': 34.0 * u.um,
            'spacing_fp': 0.738 * u.pixel,
            'width': (0.962 if use_preflightcore else 4.5) * u.pixel,
            'CDELT': [0.600737, 0.600737]*u.arcsec,
        },
    }


@u.quantity_input
@validate_channel('channel', valid_channels=[94, 131, 171, 193, 211, 304, 335]*u.angstrom)
def psf(channel: u.angstrom, use_preflightcore=False, diffraction_orders=None, use_gpu=True):
    r"""
    Calculate the composite PSF for a given channel, including diffraction and
    core effects.

    .. note:: This function has been adapted from
              `aia_calc_psf.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/psf/PRO/aia_calc_psf.pro>`_.

    .. note:: If the `cupy` package is installed
              and your machine has an NVIDIA GPU, the PSF calculation will
              automatically be accelerated with CUDA. This can lead to
              several orders of magnitude in performance increase compared to
              pure `numpy` on a CPU.

    The point spread function (PSF) can be modeled as a 2D Gaussian function
    of the radial distance :math:`r` from the center,

    .. math::

        I(r, \\theta) = I_0 \exp\left(\\frac{-r^2}{2\sigma^2}\\right)

    where,

    - :math:`I_0` : the intensity of a diffraction spike
    - :math:`r` : the radial distance from the center
    - :math:`\\theta = m\lambda/d`
    - :math:`m` : diffraction order
    - :math:`\lambda` : the wavelength of light
    - :math:`\sigma` : width of Gaussian

    The intensity of a particular diffraction spike, :math:`I_0`, is given by,

    .. math::

        I_0 = \mathrm{sinc}^2\left(\\frac{\\theta w}{\lambda}\\right)

    where,

    - :math:`w` : the width of the mesh wires
    - :math:`d` : spacing between two consecutive mesh wires

    The PSF for a given filter can then be calculated as,

    .. math::

        \mathrm{PSF} = \sum_{m=-\infty}^{+\infty}I_m(r,\\theta)

    where, in practice, one can approximate the summation by simply summing
    over a sufficiently large number of diffraction orders. In this case, we
    sum from :math:`m=--100` to :math:`m=100`.

    Finally, the composite PSF of the entrance and focal plane filters is
    given by,

    .. math::

        \mathrm{PSF}_c = \left|\mathcal{F}\left\{
                            \mathcal{F}\{\mathrm{PSF}_f\}
                            \mathcal{F}\{\mathrm{PSF}_e\}
                          \\right\}\\right|

    where :math:`\mathcal{F}` denotes the Fourier transform,
    :math:`\mathrm{PSF}_f` is the PSF of the focal plane filter, and
    :math:`\mathrm{PSF}_e` is the PSF of the entrance filter. For a more
    detailed explanation of the PSF and the above calculation, see [1]_.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of channel
    use_preflightcore : `bool`, optional
        If True, use the pre-flight values of the mesh width
    diffraction_orders : array-like, optional
        The diffraction orders to sum over. If None, the full
        range from -100 to +100 in steps of 1 will be used.
    use_gpu : `bool`, optional
        If True and `~cupy` is installed, do PSF deconvolution on the GPU
        with `~cupy`.

    Returns
    -------
    `~numpy.ndarray`
        The composite PSF of the entrance and focal plane filters.

    See Also
    --------
    filter_mesh_parameters
    deconvolve

    References
    ----------
    .. [1] `Grigis, P., Su, Y., Weber M., et al., 2012,
            AIA PSF Characterization and Deconvolution
            <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/psf/DOC/psfreport.pdf>`__
    """
    meshinfo = filter_mesh_parameters(use_preflightcore=use_preflightcore)
    meshinfo = meshinfo[channel]
    angles_entrance = meshinfo['angle_arm']
    angles_focal_plane = u.Quantity([45.0, -45.0], 'deg')
    if diffraction_orders is None:
        diffraction_orders = np.arange(-100, 101, 1)
    psf_entrance = _psf(meshinfo, angles_entrance, diffraction_orders, use_gpu=use_gpu)
    psf_focal_plane = _psf(meshinfo, angles_focal_plane, diffraction_orders,
                           focal_plane=True, use_gpu=use_gpu)
    # Composite PSF
    psf = abs(np.fft.fft2(np.fft.fft2(psf_focal_plane)
                          * np.fft.fft2(psf_entrance)))
    # Center PSF in the middle of the image
    psf = np.roll(np.roll(psf, psf.shape[1]//2, axis=1),
                  psf.shape[0]//2,
                  axis=0)
    # Normalize by total number of pixels
    psf = psf/(psf.shape[0]*psf.shape[1])
    # If using cupy, cast back to a normal numpy array
    if HAS_CUPY and use_gpu:
        psf = cupy.asnumpy(psf)
    return psf


def _psf(meshinfo, angles, diffraction_orders, focal_plane=False, use_gpu=True):
    psf = np.zeros((4096, 4096), dtype=float)
    # If cupy is available, cast to a cupy array
    if HAS_CUPY and use_gpu:
        psf = cupy.array(psf)
    Nx, Ny = psf.shape
    width_x = meshinfo['width'].value
    width_y = meshinfo['width'].value
    # x and y position grids
    x = np.outer(np.ones(Ny), np.arange(Nx) + 0.5)
    y = np.outer(np.arange(Ny) + 0.5, np.ones(Nx))
    if HAS_CUPY and use_gpu:
        x = cupy.array(x)
        y = cupy.array(y)
    area_not_mesh = 0.82  # fractional area not covered by the mesh
    spacing = meshinfo['spacing_fp'] if focal_plane else meshinfo['spacing_e']
    mesh_ratio = (meshinfo['mesh_pitch'] / meshinfo['mesh_width']).decompose().value
    spacing_x = spacing * np.cos(angles)
    spacing_y = spacing * np.sin(angles)
    for order in diffraction_orders:
        if order == 0:
            continue
        intensity = np.sinc(order / mesh_ratio)**2  # I_0
        for dx, dy in zip(spacing_x.value, spacing_y.value):
            x_centered = x - (0.5*Nx + dx*order + 0.5)
            y_centered = y - (0.5*Ny + dy*order + 0.5)
            # NOTE: this step is the bottleneck and is VERY slow on a CPU
            psf += np.exp(-width_x*x_centered*x_centered
                          - width_y*y_centered*y_centered)*intensity
    # Contribution from core
    psf_core = np.exp(-width_x*(x - 0.5*Nx - 0.5)**2
                      - width_y*(y - 0.5*Ny - 0.5)**2)
    psf_total = ((1 - area_not_mesh) * psf / psf.sum()
                 + area_not_mesh * psf_core / psf_core.sum())
    return psf_total
