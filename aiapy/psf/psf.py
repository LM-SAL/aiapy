"""
Calculate the point spread function (PSF) for the AIA telescopes.
"""

import numpy as np

import astropy.units as u

from sunpy import log
from sunpy.util.decorators import deprecated

from aiapy.psf.utils import filter_mesh_parameters
from aiapy.utils.decorators import validate_channel

try:
    import cupy

    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

__all__ = ["calculate_psf", "psf"]


@u.quantity_input
@validate_channel("channel", valid_channels=[94, 131, 171, 193, 211, 304, 335] * u.angstrom)
def calculate_psf(channel: u.angstrom, *, use_preflightcore=False, diffraction_orders=None, use_gpu=True):
    r"""
    Calculate the composite PSF for a given channel, including diffraction and
    core effects.

    .. note::

        This function has been adapted from
        `aia_calc_psf.pro <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/PRO/aia_calc_psf.pro>`_.

    .. note::

        If the `~cupy` package is installed
        and your machine has an NVIDIA GPU, the PSF calculation will
        automatically be accelerated with CUDA. This can lead to
        several orders of magnitude in performance increase compared to
        pure `numpy` on a CPU.

    The point spread function (PSF) can be modeled as a 2D Gaussian function
    of the radial distance :math:`r` from the center,

    .. math::

        I(r, \theta) = I_0 \exp\left(\frac{-r^2}{2\sigma^2}\right)

    where,

    - :math:`I_0` : the intensity of a diffraction spike
    - :math:`r` : the radial distance from the center
    - :math:`\theta = m\lambda/d`
    - :math:`m` : diffraction order
    - :math:`\lambda` : the wavelength of light
    - :math:`\sigma` : width of Gaussian

    The intensity of a particular diffraction spike, :math:`I_0`, is given by,

    .. math::

        I_0 = \mathrm{sinc}^2\left(\frac{\theta w}{\lambda}\right)

    where,

    - :math:`w` : the width of the mesh wires
    - :math:`d` : spacing between two consecutive mesh wires

    The PSF for a given filter can then be calculated as,

    .. math::

        \mathrm{PSF} = \sum_{m=-\infty}^{+\infty}I_m(r,\theta)

    where, in practice, one can approximate the summation by simply summing
    over a sufficiently large number of diffraction orders. In this case, we
    sum from :math:`m=--100` to :math:`m=100`.

    Finally, the composite PSF of the entrance and focal plane filters is
    given by,

    .. math::

        \mathrm{PSF}_c = \left|\mathcal{F}\left\{
                            \mathcal{F}\{\mathrm{PSF}_f\}
                            \mathcal{F}\{\mathrm{PSF}_e\}
                          \right\}\right|

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
            <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/DOC/psfreport.pdf>`__
    """
    meshinfo = filter_mesh_parameters(use_preflightcore=use_preflightcore)
    meshinfo = meshinfo[channel]
    angles_entrance = meshinfo["angle_arm"]
    angles_focal_plane = u.Quantity([45.0, -45.0], "deg")
    if diffraction_orders is None:
        diffraction_orders = np.arange(-100, 101, 1)
    psf_entrance = _psf(meshinfo, angles_entrance, diffraction_orders, use_gpu=use_gpu)
    psf_focal_plane = _psf(
        meshinfo,
        angles_focal_plane,
        diffraction_orders,
        focal_plane=True,
        use_gpu=use_gpu,
    )
    # Composite PSF
    psf = abs(np.fft.fft2(np.fft.fft2(psf_focal_plane) * np.fft.fft2(psf_entrance)))
    # Center PSF in the middle of the image
    psf = np.roll(np.roll(psf, psf.shape[1] // 2, axis=1), psf.shape[0] // 2, axis=0)
    # Normalize by total number of pixels
    psf = psf / (psf.shape[0] * psf.shape[1])
    # If using cupy, cast back to a normal numpy array
    if HAS_CUPY and use_gpu:
        psf = cupy.asnumpy(psf)
    return psf


def _psf(meshinfo, angles, diffraction_orders, *, focal_plane=False, use_gpu=True):
    psf = np.zeros((4096, 4096), dtype=float)
    if use_gpu and not HAS_CUPY:
        log.info("cupy not installed or working, falling back to CPU")
    # If cupy is available, cast to a cupy array
    if HAS_CUPY and use_gpu:
        psf = cupy.array(psf)
    nx, ny = psf.shape
    width_x = meshinfo["width"].value
    width_y = meshinfo["width"].value
    # x and y position grids
    x = np.outer(np.ones(ny), np.arange(nx) + 0.5)
    y = np.outer(np.arange(ny) + 0.5, np.ones(nx))
    if HAS_CUPY and use_gpu:
        x = cupy.array(x)
        y = cupy.array(y)
    area_not_mesh = 0.82  # fractional area not covered by the mesh
    spacing = meshinfo["spacing_fp"] if focal_plane else meshinfo["spacing_e"]
    mesh_ratio = (meshinfo["mesh_pitch"] / meshinfo["mesh_width"]).decompose().value
    spacing_x = spacing * np.cos(angles)
    spacing_y = spacing * np.sin(angles)
    for order in diffraction_orders:
        if order == 0:
            continue
        intensity = np.sinc(order / mesh_ratio) ** 2  # I_0
        for dx, dy in zip(spacing_x.value, spacing_y.value, strict=True):
            x_centered = x - (0.5 * nx + dx * order + 0.5)
            y_centered = y - (0.5 * ny + dy * order + 0.5)
            # NOTE: this step is the bottleneck and is VERY slow on a CPU
            psf += np.exp(-width_x * x_centered * x_centered - width_y * y_centered * y_centered) * intensity
    # Contribution from core
    psf_core = np.exp(-width_x * (x - 0.5 * nx - 0.5) ** 2 - width_y * (y - 0.5 * ny - 0.5) ** 2)
    return (1 - area_not_mesh) * psf / psf.sum() + area_not_mesh * psf_core / psf_core.sum()


@deprecated("0.11.0", alternative="calculate_psf")
@u.quantity_input
def psf(channel: u.angstrom, *, use_preflightcore=False, diffraction_orders=None, use_gpu=True):
    """
    This function is deprecated.

    Use `aiapy.psf.calculate_psf` instead.
    """
    return calculate_psf(
        channel,
        use_preflightcore=use_preflightcore,
        diffraction_orders=diffraction_orders,
        use_gpu=use_gpu,
    )
