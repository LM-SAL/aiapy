import pytest

import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from aiapy.utils import check_quality_flag, sdo_location
from aiapy.utils.utils import _QUALITY_FLAG_MESSAGES


@pytest.mark.remote_data
def test_sdo_location(aia_171_test_map) -> None:
    # Confirm that the queried location matches AIAMap's interpretation of the FITS file
    result = sdo_location(aia_171_test_map.reference_date)
    observer = aia_171_test_map.observer_coordinate.transform_to(result)
    assert_quantity_allclose(result.cartesian.xyz, observer.cartesian.xyz, rtol=0, atol=1 * u.cm)


@pytest.mark.parametrize(
    "bits",
    [
        [],  # Nominal
        [16],  # Single message
        [12, 13, 14, 17, 21],  # Multiple messages
        [4, 26],  # Empty bits
    ],
)
def test_check_quality_flag(bits):
    quality = 0
    for b in bits:
        quality = quality | (1 << b)
    messages = ["nominal"]
    if bits:
        messages = [_QUALITY_FLAG_MESSAGES.get(b, "(empty)") for b in bits]
    assert messages == check_quality_flag(quality)


def test_check_quality_flag_hand_set_bits():
    # Set by hand at JSOC, so missing from the SDO User Guide table
    assert "(empty)" not in check_quality_flag(sum(1 << b for b in (5, 6, 7, 24, 25)))
