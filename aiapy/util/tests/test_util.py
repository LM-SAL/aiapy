import pytest

from astropy.tests.helper import assert_quantity_allclose

import aiapy.util


@pytest.mark.remote_data()
def test_sdo_location(aia_171_map):
    # Confirm that the queried location matches AIAMap's interpretation of the FITS file
    result = aiapy.util.sdo_location(aia_171_map.date)
    aia_171_map.observer_coordinate.transform_to(result)
    assert_quantity_allclose(result.cartesian.xyz, result.cartesian.xyz)


@pytest.mark.remote_data()
def test_sdo_location_raises_error():
    # Confirm that an error is raised for a time without records
    with pytest.raises(ValueError, match="No DRMS records near this time"):
        aiapy.util.sdo_location("2001-01-01")
