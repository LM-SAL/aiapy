"""
This package contains functions for deconvolving AIA images with the instrument
point spread function (PSF).
"""

try:
    from functools import partial

    import jax

    jax.config.update("jax_enable_x64", True)  # NOQA: FBT003

    import jax.numpy as np
    from jax import jit, lax

    _jit_iterations = partial(jit, static_argnames=("iterations",))

except ImportError:
    import numpy as np  # NOQA: F401

    class _LaxShim:
        @staticmethod
        def fori_loop(lower, upper, body_fun, init_val):
            acc = init_val
            for i in range(lower, upper):
                acc = body_fun(i, acc)
            return acc

        @staticmethod
        def dynamic_index_in_dim(arr, i, keepdims=False):  # NOQA: FBT002 ARG004
            return arr[i]

    lax = _LaxShim()

    def _jit_iterations(fun=None, **kwargs):  # NOQA: ARG001
        return fun


from .deconvolve import *
from .psf import *
from .utils import *
