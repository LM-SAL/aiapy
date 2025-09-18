import numpy as np


def test_psf(psf) -> None:
    assert isinstance(psf, np.ndarray)
    assert psf.shape == (4096, 4096)
