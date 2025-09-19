"""
Calculate the point spread function (PSF) for the AIA telescopes.
"""

from numpy import asarray

import astropy.units as u

from sunpy.util.decorators import deprecated

from aiapy.psf import np
from aiapy.psf.utils import filter_mesh_parameters
from aiapy.utils.decorators import validate_channel

__all__ = ["calculate_psf", "psf"]


@u.quantity_input
@validate_channel("channel", valid_channels=[94, 131, 171, 193, 211, 304, 335] * u.angstrom)
def calculate_psf(channel: u.angstrom, *, use_preflightcore=False, diffraction_orders=None):
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
    mesh_info = filter_mesh_parameters(use_preflightcore=use_preflightcore)[channel]
    diffraction_orders = np.arange(-100, 101, 1) if diffraction_orders is None else np.asarray(diffraction_orders)
    psf_e = _psf(mesh_info, mesh_info["angle_arm"], diffraction_orders, focal_plane=False)
    psf_f = _psf(mesh_info, u.Quantity([45.0, -45.0], "deg"), diffraction_orders, focal_plane=True)
    # Composite PSF
    psf = np.abs(np.fft.fft2(np.fft.fft2(psf_e) * np.fft.fft2(psf_f)))
    # Center PSF at pixel (0,0)
    psf = np.fft.fftshift(psf)
    # Normalize by total number of pixels and always return a numpy array
    return asarray(psf / (psf.shape[0] * psf.shape[1]))


def _psf(mesh_info, angles, diffraction_orders, *, focal_plane=False):
    ny, nx = 4096, 4096
    width = mesh_info["width"].to_value("pixel")
    spacing = (mesh_info["spacing_fp"] if focal_plane else mesh_info["spacing_e"]).to_value("pixel")
    mesh_ratio = (mesh_info["mesh_pitch"] / mesh_info["mesh_width"]).to_value()
    area_not_mesh = 0.82
    x = np.arange(nx) + 0.5
    y = np.arange(ny) + 0.5
    cx0 = 0.5 * nx + 0.5
    cy0 = 0.5 * ny + 0.5
    ang = np.asarray(angles.to_value(u.rad))
    dxs = spacing * np.cos(ang)
    dys = spacing * np.sin(ang)
    psf_spikes = np.zeros((ny, nx))
    for order in diffraction_orders:
        if order == 0:
            continue
        i0 = np.sinc(order / mesh_ratio) ** 2
        for dx, dy in zip(dxs, dys, strict=True):
            gx = np.exp(-width * (x - (cx0 + dx * order)) ** 2)
            gy = np.exp(-width * (y - (cy0 + dy * order)) ** 2)
            psf_spikes = psf_spikes + i0 * (gy[:, None] * gx[None, :])
    gx0 = np.exp(-width * (x - cx0) ** 2)
    gy0 = np.exp(-width * (y - cy0) ** 2)
    psf_core = gy0[:, None] * gx0[None, :]
    psf_spikes = psf_spikes / psf_spikes.sum()
    psf_core = psf_core / psf_core.sum()
    return (1.0 - area_not_mesh) * psf_spikes + area_not_mesh * psf_core


@deprecated("0.11.0", alternative="calculate_psf")
@u.quantity_input
def psf(channel: u.angstrom, *, use_preflightcore=False, diffraction_orders=None):
    """
    This function is deprecated.

    Use `aiapy.psf.calculate_psf` instead.
    """
    return calculate_psf(
        channel,
        use_preflightcore=use_preflightcore,
        diffraction_orders=diffraction_orders,
    )
