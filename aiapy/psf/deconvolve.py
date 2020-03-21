"""
Deconvolve an AIA image with the channel point spread function
"""
import copy

import numpy as np
try:
    import cupy
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False

from .psf import psf as calculate_psf

__all__ = ['deconvolve']


def deconvolve(smap, psf=None, iterations=25):
    """
    Deconvolve an AIA image with the point spread function

    Perform image deconvolution on an AIA image with the instrument
    point spread function using the Richardson-Lucy deconvolution
    algorithm [1]_.

    .. note:: If the `cupy` package is installed
              and your machine has an NVIDIA GPU, the deconvolution will
              automatically be accelerated with CUDA. This can lead to more
              than an order of magnitude in performance increase compared to
              pure `numpy` on a CPU. For more information on PSF deconvolution
              on a GPU, see [2]_.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        An AIA image
    psf : `~numpy.ndarray`, optional
        The point spread function. If None, it will be calculated
    iterations: `int`
        Number of iterations in the Richardson-Lucy algorithm

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
    .. [2] Cheung, M., 2015, *GPU Technology Conference Silicon Valley*, `GPU-Accelerated Image Processing for NASA's Solar Dynamics Observatory <https://on-demand-gtc.gputechconf.com/gtcnew/sessionview.php?sessionName=s5209-gpu-accelerated+imaging+processing+for+nasa%27s+solar+dynamics+observatory>`_
    """
    # TODO: do we need a check to make sure this is a full-frame image?
    img = smap.data
    if psf is None:
        psf = calculate_psf(smap.wavelength)
    if HAS_CUPY:
        img = cupy.array(img)
        psf = cupy.array(psf)

    # Center PSF at pixel (0,0)
    psf = np.roll(np.roll(psf, psf.shape[0]//2, axis=0),
                  psf.shape[1]//2,
                  axis=1)
    # Convolution requires FFT of the PSF
    psf = np.fft.rfft2(psf)
    psf_conj = psf.conj()

    img_decon = np.copy(img)
    for _ in range(iterations):
        ratio = img/np.fft.irfft2(np.fft.rfft2(img_decon)*psf)
        img_decon = img_decon*np.fft.irfft2(np.fft.rfft2(ratio)*psf_conj)

    return smap._new_instance(
        cupy.asnumpy(img_decon) if HAS_CUPY else img_decon,
        copy.deepcopy(smap.meta),
        plot_settings=copy.deepcopy(smap.plot_settings),
        mask=smap.mask,
    )
