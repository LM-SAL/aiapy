from contextlib import nullcontext

import numpy as np
import pytest

import astropy.units as u

from aiapy.calibrate import estimate_error
from aiapy.calibrate.util import get_error_table
from aiapy.conftest import CHANNELS, RANDOM_GENERATOR
from aiapy.tests.data import get_test_filepath

table_local = get_error_table(get_test_filepath("aia_V3_error_table.txt"))


@pytest.mark.parametrize("channel", CHANNELS)
def test_error_all_channels(channel) -> None:
    intensity = 10.0 * u.DN / u.pix
    error = estimate_error(intensity, channel, error_table=table_local)
    assert error.unit == intensity.unit


@pytest.mark.parametrize(
    "counts",
    [
        1,
        RANDOM_GENERATOR.standard_normal(1),
        RANDOM_GENERATOR.standard_normal(10),
        RANDOM_GENERATOR.standard_normal((100, 200)),
        RANDOM_GENERATOR.standard_normal((10, 10, 5)),
    ],
)
def test_counts_shapes(counts) -> None:
    counts = np.abs(counts) * 1000 * u.DN / u.pix
    errors = estimate_error(counts, 171 * u.angstrom, error_table=table_local)
    if counts.shape == ():
        assert errors.shape == (1,)
    else:
        assert counts.shape == errors.shape


@pytest.mark.parametrize(
    ("include_preflight", "include_eve", "include_chianti", "expectation"),
    [
        (False, True, False, nullcontext()),
        (True, False, False, nullcontext()),
        (False, False, False, nullcontext()),
        (False, False, True, nullcontext()),
        (
            True,
            True,
            False,
            pytest.raises(ValueError, match="Cannot include both EVE and pre-flight correction."),
        ),
    ],
)
def test_flags(include_preflight, include_eve, include_chianti, expectation) -> None:
    with expectation:
        errors = estimate_error(
            1 * u.DN / u.pix,
            94 * u.angstrom,
            error_table=table_local,
            include_eve=include_eve,
            include_preflight=include_preflight,
            include_chianti=include_chianti,
        )
        assert isinstance(errors, u.Quantity)
