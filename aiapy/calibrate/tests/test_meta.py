"""
Tests for function that update FITS keywords
"""
import pytest

import sunpy.data.test
from sunpy.map import Map

from aiapy.calibrate import fix_observer_location, update_pointing


@pytest.fixture
def smap():
    return Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits'))


def test_fix_observer_location(smap):
    smap_fixed = fix_observer_location(smap)
    # NOTE: AIAMap already fixes the .observer_coordinate property with HAE
    assert smap_fixed.meta['hgln_obs'] == smap_fixed.observer_coordinate.lon.value
    assert smap_fixed.meta['hglt_obs'] == smap_fixed.observer_coordinate.lat.value
    assert smap_fixed.meta['dsun_obs'] == smap_fixed.observer_coordinate.radius.value


@pytest.mark.remote_data
def test_fix_pointing(smap):
    keys = ['CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CROTA2']
    # Remove keys to at least test that they get set
    for k in keys:
        _ = smap.meta.pop(k)

    aia_map_updated = update_pointing(smap)
    # FIXME: how do we check these values are accurate?
    assert all([k in aia_map_updated.meta for k in keys])
