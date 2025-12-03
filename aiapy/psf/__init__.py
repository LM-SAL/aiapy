"""
This package contains functions for deconvolving AIA images with the instrument
point spread function (PSF).
"""

from sunpy import log

# The goal here is to use jax if it is installed, otherwise fall back to numpy.
# If JAX is not installed we need to provide fake versions of various jax
# functions and objects that we use. These are basically no-ops or simple
# implementations using numpy.
try:
    from functools import partial

    import jax

    # If we do not do this, jax will use 32-bit floats by default which
    # results in significant errors for our calculations.
    jax.config.update("jax_enable_x64", True)  # NOQA: FBT003

    import jax.numpy as np
    from jax import jit, lax

    _jit_over_iterations = partial(jit, static_argnames=("iterations",))
    log.debug("Using jax for the PSF module.")
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

    # From the JAX documentation:
    # jax.lax is a library of primitives operations that underpins libraries such as jax.numpy.
    # Transformation rules, such as JVP and batching rules, are typically defined as
    # transformations on jax.lax primitives.
    #
    # These are low-level functions that are not typically called directly by users.
    # However, we need lax.fori_loop and lax.dynamic_index_in_dim as they are used to
    # implement loops and indexing in jax.jit compiled functions.
    lax = _LaxShim()

    # We need to provide a no-op version of the jit decorator when jax is not
    # installed.
    def jit(func):
        return func

    def _jit_over_iterations(fun=None, **kwargs):  # NOQA: ARG001
        return fun

    log.debug("Using numpy for the PSF module.")

from .deconvolve import *
from .psf import *
from .utils import *

__all__ = ["calculate_psf", "deconvolve", "filter_mesh_parameters", "psf"]  # NOQA: F405
