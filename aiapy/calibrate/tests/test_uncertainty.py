from contextlib import nullcontext

import astropy.units as u
import numpy as np
import pytest
from numpy.random import default_rng

from aiapy.calibrate import estimate_error
from aiapy.calibrate.util import get_error_table
from aiapy.tests.data import get_test_filepath

# These are not fixtures so that they can be easily used in the parametrize mark
RANDOM_GENERATOR = default_rng()
CHANNELS = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500] * u.angstrom
table_local = get_error_table(get_test_filepath("aia_V3_error_table.txt"))


@pytest.mark.parametrize("channel", CHANNELS)
def test_error_all_channels(channel):
    intensity = 10.0 * u.ct / u.pix
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
def test_counts_shapes(counts):
    counts = np.abs(counts) * 1000 * u.ct / u.pix
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
def test_flags(include_preflight, include_eve, include_chianti, expectation):
    with expectation:
        errors = estimate_error(
            1 * u.ct / u.pix,
            94 * u.angstrom,
            error_table=table_local,
            include_eve=include_eve,
            include_preflight=include_preflight,
            include_chianti=include_chianti,
        )
        assert isinstance(errors, u.Quantity)
