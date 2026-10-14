import numpy as np

from aiapy.utils.utils import detector_dimensions


def test_psf(psf) -> None:
    # We always return a numpy array, even if JAX is used as a backend
    assert isinstance(psf, np.ndarray)
    np.testing.assert_equal(psf.shape, detector_dimensions().value)
    # Normalized, non-negative and centered on the detector
    np.testing.assert_allclose(psf.sum(), 1, rtol=1e-10)
    assert psf.min() >= 0
    assert np.unravel_index(psf.argmax(), psf.shape) == (2048, 2048)
