import numpy as np


def test_smoke_psf(psf) -> None:
    # We always return a numpy array, even if JAX is used as a backend
    assert isinstance(psf, np.ndarray)
    assert psf.shape == (4096, 4096)
