"""
Deconvolve an AIA image with the channel point spread function.
"""

import copy
import warnings

import numpy as np

try:
    import jax.numpy as jnp
    from jax import jit as jax_jit
except ImportError:
    jnp = np

    def jax_jit(func):
        return func


from aiapy.util import AIApyUserWarning
from .psf import psf as calculate_psf

__all__ = ["deconvolve"]


@jax_jit
def _deconvolve_inner(img, psf, psf_conj, iterations):
    img_decon = img.copy()
    for _ in range(iterations):
        ratio = img / jnp.fft.irfft2(jnp.fft.rfft2(img_decon) * psf)
        img_decon = img_decon * jnp.fft.irfft2(jnp.fft.rfft2(ratio) * psf_conj)
    return img_decon


def deconvolve(smap, *, psf=None, iterations=25, clip_negative=True):
    """
    Deconvolve an AIA image with the point spread function.

    Perform image deconvolution on an AIA image with the instrument
    point spread function using the Richardson-Lucy deconvolution
    algorithm [1]_.

    .. note::

        If the `~jax` package is installed it will be used to accelerate the computation.
        `~jax` can use CPUs or GPUs, see their documentation for instructions.
        For more information on PSF deconvolution on a GPU, see [2]_.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        An AIA image
    psf : `~numpy.ndarray`, optional
        The point spread function. If None, it will be calculated
    iterations : `int`, optional
        Number of iterations in the Richardson-Lucy algorithm
    clip_negative : `bool`, optional
        If the image has negative intensity values, set them to zero.

    Returns
    -------
    `~sunpy.map.Map`
        Deconvolved AIA image

    See Also
    --------
    psf

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Richardson%E2%80%93Lucy_deconvolution
    .. [2] Cheung, M., 2015, *GPU Technology Conference Silicon Valley*, `GPU-Accelerated Image Processing for NASA's Solar Dynamics Observatory <https://on-demand.gputechconf.com/gtc/2015/presentation/S5209-Mark-Cheung.pdf>`__
    """
    # TODO: do we need a check to make sure this is a full-frame image?
    img = smap.data
    if clip_negative:
        img = jnp.where(img < 0, 0, img)
    if jnp.any(img < 0):
        warnings.warn(
            "Image contains negative intensity values. Consider setting clip_negative to True",
            AIApyUserWarning,
            stacklevel=3,
        )
    if psf is None:
        psf = calculate_psf(smap.wavelength)
    # Center PSF at pixel (0,0)
    psf = jnp.roll(jnp.roll(psf, psf.shape[0] // 2, axis=0), psf.shape[1] // 2, axis=1)
    # Convolution requires FFT of the PSF
    psf = jnp.fft.rfft2(psf)
    psf_conj = psf.conj()
    img_decon = _deconvolve_inner(img, psf, psf_conj, iterations)
    return smap._new_instance(
        img_decon,
        copy.deepcopy(smap.meta),
        plot_settings=copy.deepcopy(smap.plot_settings),
        mask=smap.mask,
    )
