"""
Deconvolve an AIA image with the channel point spread function.
"""

import copy
import warnings

from aiapy.psf import _jit_iterations, lax, np
from aiapy.psf.psf import calculate_psf
from aiapy.utils import AIApyUserWarning

__all__ = ["deconvolve"]


@_jit_iterations
def _rl_deconvolve(img, psf_fft, psf_conj, *, iterations: int):
    """
    Richardson-Lucy deconvolution core.
    """

    def body(_, current):
        est = np.fft.irfft2(np.fft.rfft2(current) * psf_fft)
        ratio = img / est
        update = np.fft.irfft2(np.fft.rfft2(ratio) * psf_conj)
        return current * update

    return lax.fori_loop(0, iterations, body, img)


def deconvolve(
    smap,
    *,
    psf: np.ndarray | None = None,
    iterations: int = 25,
    clip_negative: bool = True,
):
    """
    Deconvolve an AIA image with the point spread function.

    Performs Richardson-Lucy deconvolution.

    .. note::

        If the jax package is installed it will be used to accelerate the computation.
        jax can use CPUs or GPUs, `see their documentation for instructions <https://docs.jax.dev/en/latest/installation.html>`__.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        An AIA image.
    psf : array-like, optional
        The point spread function.
        Defaults to `None` and it will be calculated with `aiapy.psf.calculate_psf`.
    iterations : `int`, optional
        Number of Richardson-Lucy iterations, defaults to 25.
    clip_negative : `bool`, optional
        If the image has negative intensity values, set them to zero.
        Defaults to `True`.

    Returns
    -------
    `~sunpy.map.Map`
        Deconvolved AIA image

    See Also
    --------
    calculate_psf

    References
    ----------
    https://en.wikipedia.org/wiki/Richardson%E2%80%93Lucy_deconvolution
    Cheung, M., 2015, *GPU Technology Conference Silicon Valley*, `GPU-Accelerated Image Processing for NASA's Solar Dynamics Observatory <https://on-demand.gputechconf.com/gtc/2015/presentation/S5209-Mark-Cheung.pdf>`__
    """
    # TODO: Should check to make sure this is a full-frame image?
    img = smap.data.astype(np.float64)
    if clip_negative:
        img = np.where(img < 0, 0, img)
    if np.any(img < 0):
        warnings.warn(
            "Image contains negative intensity values. Consider setting clip_negative to True",
            AIApyUserWarning,
            stacklevel=2,
        )
    if psf is None:
        psf = calculate_psf(smap.wavelength)
    # Center PSF at pixel (0,0)
    psf = np.roll(np.roll(psf, psf.shape[0] // 2, axis=0), psf.shape[1] // 2, axis=1)
    # Convolution requires FFT of the PSF
    psf = np.fft.rfft2(psf)
    psf_conj = np.conj(psf)
    img_decon = _rl_deconvolve(img, psf, psf_conj, iterations=iterations)
    return smap._new_instance(
        img_decon,
        copy.deepcopy(smap.meta),
        plot_settings=copy.deepcopy(smap.plot_settings),
        mask=smap.mask,
    )
