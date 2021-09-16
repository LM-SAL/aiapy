import numpy as np
import pytest

import sunpy.data.test
import sunpy.map

import aiapy.psf
from aiapy.util import AiapyUserWarning


def test_deconvolve(aia_171_map):
    # Skip this test if cupy is not installed because it is too
    # slow. This is mostly for the benefit of the CI.
    try:
        import cupy  # NOQA
    except ImportError:
        pytest.skip('Cannot import cupy. Skipping deconvolution test with full PSF')
    map_decon = aiapy.psf.deconvolve(aia_171_map)
    assert isinstance(map_decon, sunpy.map.GenericMap)
    assert map_decon.data.shape == aia_171_map.data.shape


def test_deconvolve_specify_psf(aia_171_map, psf):
    map_decon = aiapy.psf.deconvolve(aia_171_map, psf=psf, iterations=1)
    assert isinstance(map_decon, sunpy.map.GenericMap)
    assert map_decon.data.shape == aia_171_map.data.shape


def test_deconvolve_negative_pixels(aia_171_map, psf):
    aia_171_map_neg = aia_171_map._new_instance(
        np.where(aia_171_map.data < 1, -1, aia_171_map.data),
        aia_171_map.meta,
    )
    with pytest.warns(AiapyUserWarning):
        _ = aiapy.psf.deconvolve(
            aia_171_map_neg,
            psf=psf,
            iterations=1,
            clip_negative=False,
        )
