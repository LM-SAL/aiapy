from contextlib import nullcontext

import numpy as np
import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from aiapy.calibrate import estimate_error
from aiapy.calibrate.utils import get_error_table
from aiapy.tests.data import get_test_filepath

table_local = get_error_table(get_test_filepath("aia_V3_error_table.txt"))


@pytest.mark.parametrize(
    ("channel", "counts", "include_eve", "include_chianti", "expected"),
    # Reference values from SSW: aia_bp_estimate_error(counts, channel, n_sample=1),
    # see test_idl.py::test_error_consistent
    [
        (94, 10, False, False, 4.603076840668012),
        (131, 10, False, False, 4.020984506447918),
        (171, 10, False, False, 3.5578527795700317),
        (193, 10, False, False, 3.3254046973404905),
        (211, 10, False, False, 3.1997369268573),
        (304, 10, False, False, 2.7474199535371726),
        (335, 10, False, False, 2.6829678345247077),
        (171, 1000, True, False, 38.598564688738406),
        (171, 1000, False, True, 252.24759502526624),
    ],
)
def test_error_matches_idl(channel, counts, include_eve, include_chianti, expected) -> None:
    error = estimate_error(
        counts * u.DN / u.pix,
        channel * u.angstrom,
        error_table=table_local,
        include_eve=include_eve,
        include_chianti=include_chianti,
        compare_idl=True,
    )
    assert_quantity_allclose(error, [expected] * u.DN / u.pix, rtol=1e-4)


@pytest.mark.parametrize("shape", [(), (1,), (10,), (100, 200), (10, 10, 5)])
def test_counts_shapes(shape) -> None:
    counts = np.full(shape, 1000.0) * u.DN / u.pix
    errors = estimate_error(counts, 171 * u.angstrom, error_table=table_local)
    assert errors.shape == (counts.shape or (1,))


@pytest.mark.parametrize(
    ("include_preflight", "include_eve", "expectation"),
    [
        (True, False, nullcontext()),
        (True, True, pytest.raises(ValueError, match=r"Cannot include both EVE and pre-flight correction.")),
    ],
)
def test_flags(include_preflight, include_eve, expectation) -> None:
    with expectation:
        errors = estimate_error(
            1 * u.DN / u.pix,
            94 * u.angstrom,
            error_table=table_local,
            include_eve=include_eve,
            include_preflight=include_preflight,
        )
        assert isinstance(errors, u.Quantity)
