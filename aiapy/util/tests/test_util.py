"""
Test utilities
"""
from astropy.tests.helper import assert_quantity_allclose
import pytest
import sunpy.data.test
from sunpy.map import Map

import aiapy.util


@pytest.fixture
def smap():
    return Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits'))


@pytest.mark.remote_data
def test_sdo_location(smap):
    # Confirm that the queried location matches AIAMap's interpretation of the FITS file
    result = aiapy.util.sdo_location(smap.date)
    check = smap.observer_coordinate.transform_to(result)

    assert_quantity_allclose(result.cartesian.xyz, result.cartesian.xyz)


@pytest.mark.remote_data
def test_sdo_location():
    # Confirm that an error is raised for a time without records
    with pytest.raises(ValueError):
        aiapy.util.sdo_location('2001-01-01')
