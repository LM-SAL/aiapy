from itertools import permutations

import numpy as np
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from aiapy.analysis import isolate_fe18_intensity
from aiapy.tests.helpers import make_aia_map


@pytest.mark.parametrize("order", list(permutations(range(3))))
def test_isolate_fe18_intensity(order) -> None:
    """Check the published coefficients for every input order."""
    # (I_171, I_193, I_94warm), in DN/s. Single-channel inputs test each weight.
    intensities = np.array(
        [
            (0.0, 0.0, -0.028509),
            (100.0, 0.0, 0.0753497803),
            (0.0, 100.0, 0.2099309435),
            (1000.0, 3000.0, 14.4060611497),
        ],
    )
    i_171 = intensities[:, 0].reshape(2, 2)
    i_193 = intensities[:, 1].reshape(2, 2)
    expected = intensities[:, 2].reshape(2, 2)
    i_94 = np.full((2, 2), 10.0)
    maps = (
        make_aia_map(i_94, 94),
        make_aia_map(i_171, 171),
        make_aia_map(i_193, 193),
    )
    result = isolate_fe18_intensity(*(maps[index] for index in order))
    assert_quantity_allclose(
        maps[0].quantity - result.quantity,
        expected * u.DN / u.s,
        rtol=1e-6,
    )
    assert result.meta == maps[0].meta
    assert result.mask is None
    for smap, before in zip(maps, (i_94, i_171, i_193), strict=True):
        np.testing.assert_array_equal(smap.data, before)


def test_composite_intensity_is_limited() -> None:
    """Cap the scaled composite at 30."""
    intensity = 116.54 * np.array([[20.0, 30.0], [40.0, 1000.0]])
    result = isolate_fe18_intensity(
        make_aia_map(np.full((2, 2), 100.0), 94),
        make_aia_map(intensity, 171),
        make_aia_map(intensity, 193),
    )
    assert_quantity_allclose(
        (100.0 - result.data) * u.DN / u.s,
        np.array([[14.159691, 16.222791], [16.222791, 16.222791]]) * u.DN / u.s,
        rtol=1e-6,
    )


@pytest.mark.parametrize("masked_indices", [(0,), (1,), (2,), (0, 1, 2)])
def test_output_combines_input_masks(masked_indices) -> None:
    maps = [make_aia_map(np.full((2, 2), 1.0), w) for w in (94, 171, 193)]
    expected = np.zeros((2, 2), dtype=bool)
    for index in masked_indices:
        mask = np.zeros((2, 2), dtype=bool)
        mask.flat[index] = True
        maps[index].mask = mask
        expected |= mask
    masks_before = [None if smap.mask is None else smap.mask.copy() for smap in maps]

    result = isolate_fe18_intensity(*maps)

    np.testing.assert_array_equal(result.mask, expected)
    for smap, mask_before in zip(maps, masks_before, strict=True):
        np.testing.assert_array_equal(smap.mask, mask_before)


@pytest.mark.parametrize("scalar_mask", [False, True])
def test_output_combines_scalar_and_array_masks(scalar_mask) -> None:
    maps = [make_aia_map(np.full((2, 2), 1.0), w) for w in (94, 171, 193)]
    maps[0].mask = scalar_mask
    maps[1].mask = np.array([[True, False], [False, False]])
    result = isolate_fe18_intensity(*maps)
    expected = np.full((2, 2), scalar_mask)
    expected[0, 0] = True
    np.testing.assert_array_equal(result.mask, expected)


@pytest.mark.parametrize("index", [0, 1, 2])
def test_inputs_must_be_normalized_by_exposure_time(index) -> None:
    wavelengths = (94, 171, 193)
    maps = [make_aia_map(np.full((8, 8), 1.0), w) for w in wavelengths]
    maps[index] = make_aia_map(np.full((8, 8), 1.0), wavelengths[index], unit=u.DN)
    with pytest.raises(ValueError, match="must have units of DN/s"):
        isolate_fe18_intensity(*maps)


def test_inputs_must_have_same_dimensions() -> None:
    m_94 = make_aia_map(np.full((8, 8), 1.0), 94)
    m_171 = make_aia_map(np.full((4, 4), 1.0), 171)
    m_193 = make_aia_map(np.full((8, 8), 1.0), 193)
    with pytest.raises(ValueError, match="must have the same dimensions"):
        isolate_fe18_intensity(m_94, m_171, m_193)


@pytest.mark.parametrize("wavelengths", [(94, 171, 211), (94, 171, 171)])
def test_inputs_must_include_each_wavelength_once(wavelengths) -> None:
    maps = [make_aia_map(np.ones((2, 2)), wavelength) for wavelength in wavelengths]
    with pytest.raises(ValueError, match="one map each at 94, 171, and 193 Å"):
        isolate_fe18_intensity(*maps)
