import jax.numpy as jnp
import numpy as np


def test_smoke_psf(psf) -> None:
    assert isinstance(psf, jnp.ndarray | np.ndarray)
    assert psf.shape == (4096, 4096)
