import pytest

from astropy.tests.helper import assert_quantity_allclose

from aiapy.utils.utils import _QUALITY_FLAG_MESSAGES, check_quality_flag, sdo_location


@pytest.mark.remote_data
def test_sdo_location(aia_171_map) -> None:
    # Confirm that the queried location matches AIAMap's interpretation of the FITS file
    result = sdo_location(aia_171_map.date)
    aia_171_map.observer_coordinate.transform_to(result)
    assert_quantity_allclose(result.cartesian.xyz, result.cartesian.xyz)


@pytest.mark.remote_data
def test_sdo_location_raises_error() -> None:
    # Confirm that an error is raised for a time without records
    with pytest.raises(RuntimeError, match="No data found for this query"):
        sdo_location("2001-01-01")


@pytest.mark.parametrize(
    "bits",
    [
        [],  # Nominal
        [16],  # Single message
        [12, 13, 14, 17, 21],  # Multiple messages
        [4, 5],  # Empty bits
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
