import numpy as np
import pytest

import sunpy.map

import aiapy.psf
from aiapy.utils import AIApyUserWarning


def test_deconvolve(aia_171_map) -> None:
    pytest.importorskip(modname="jax", reason="Cannot import jax. Skipping deconvolution test with full PSF")
    map_decon = aiapy.psf.deconvolve(aia_171_map)
    assert isinstance(map_decon, sunpy.map.GenericMap)
    assert map_decon.data.shape == aia_171_map.data.shape


def test_deconvolve_delta_psf(aia_171_positive_map) -> None:
    # A delta-function PSF centered like calculate_psf leaves the image unchanged
    psf = np.zeros(aia_171_positive_map.data.shape)
    psf[64, 64] = 1
    map_decon = aiapy.psf.deconvolve(aia_171_positive_map, psf=psf, iterations=5)
    np.testing.assert_allclose(map_decon.data, aia_171_positive_map.data, rtol=1e-10)


def test_deconvolve_matches_idl(aia_171_positive_map, psf_128) -> None:
    # Reference values from SSW: aia_deconvolve_richardsonlucy(image, psf, niter=25)
    # See test_idl.py::test_deconvolve_consistent for the full comparison.
    map_decon = aiapy.psf.deconvolve(aia_171_positive_map, psf=psf_128, iterations=25)
    np.testing.assert_allclose(
        [map_decon.data[50, 70], map_decon.data[64, 64], map_decon.data[100, 40]],
        [11721.527399200102, 130.73597027588772, 365.912288319648],
        rtol=1e-10,
    )


@pytest.mark.parametrize("value", [0, -1])
def test_deconvolve_non_positive_image(aia_171_test_map, psf_128, value) -> None:
    # After clipping the image and every estimate are zero, which must not become NaN
    smap = aia_171_test_map._new_instance(np.full((128, 128), value, dtype=float), aia_171_test_map.meta)
    np.testing.assert_array_equal(aiapy.psf.deconvolve(smap, psf=psf_128, iterations=2).data, 0)


def test_deconvolve_negative_pixels(aia_171_test_map, psf_128) -> None:
    negative = aia_171_test_map.data < 0
    assert negative.any()
    # Clipped pixels start at zero and the multiplicative update keeps them there
    map_decon = aiapy.psf.deconvolve(aia_171_test_map, psf=psf_128, iterations=5)
    np.testing.assert_array_equal(map_decon.data[negative], 0)
    with pytest.warns(AIApyUserWarning, match="Image contains negative intensity values."):
        aiapy.psf.deconvolve(aia_171_test_map, psf=psf_128, iterations=1, clip_negative=False)
