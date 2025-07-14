import pytest

from astropy.tests.helper import assert_quantity_allclose

import aiapy.util
from aiapy.util.util import _QUALITY_FLAG_MESSAGES


@pytest.mark.remote_data
def test_sdo_location(aia_171_map) -> None:
    # Confirm that the queried location matches AIAMap's interpretation of the FITS file
    result = aiapy.util.sdo_location(aia_171_map.date)
    aia_171_map.observer_coordinate.transform_to(result)
    assert_quantity_allclose(result.cartesian.xyz, result.cartesian.xyz)


@pytest.mark.remote_data
def test_sdo_location_raises_error() -> None:
    # Confirm that an error is raised for a time without records
    with pytest.raises(RuntimeError, match="No data found for this query"):
        aiapy.util.sdo_location("2001-01-01")


@pytest.mark.parametrize(
    "bits",
    [
        [],  # nominal
        [16],  # single message
        [12, 13, 14, 17, 21],  # multiple messages
        [4, 5],  # empty bits
    ],
)
def test_check_quality_flag(bits):
    quality = 0
    for b in bits:
        quality = quality | (1 << b)
    messages = ["nominal"]
    if bits:
        messages = [_QUALITY_FLAG_MESSAGES.get(b, "(empty)") for b in bits]
    assert messages == aiapy.util.check_quality_flag(quality)
