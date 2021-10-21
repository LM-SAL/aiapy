import pytest

import astropy.units as u

from aiapy.calibrate import fix_observer_location, update_pointing
from aiapy.calibrate.util import get_pointing_table
from aiapy.util.exceptions import AiapyUserWarning


def test_fix_observer_location(aia_171_map):
    smap_fixed = fix_observer_location(aia_171_map)
    # NOTE: AIAMap already fixes the .observer_coordinate property with HAE
    assert smap_fixed.meta['hgln_obs'] == smap_fixed.observer_coordinate.lon.value
    assert smap_fixed.meta['hglt_obs'] == smap_fixed.observer_coordinate.lat.value
    assert smap_fixed.meta['dsun_obs'] == smap_fixed.observer_coordinate.radius.value


@pytest.mark.remote_data
def test_fix_pointing(aia_171_map):
    keys = ['CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CROTA2']
    # Remove keys to at least test that they get set
    for k in keys:
        _ = aia_171_map.meta.pop(k)
    aia_map_updated = update_pointing(aia_171_map)
    # FIXME: how do we check these values are accurate?
    assert all([k in aia_map_updated.meta for k in keys])
    # Check the case where we have specified the pointing
    # table ahead of time
    ptable = get_pointing_table(aia_171_map.date-6*u.h,
                                aia_171_map.date+6*u.h)
    aia_map_updated2 = update_pointing(aia_171_map, pointing_table=ptable)
    for k in keys:
        assert aia_map_updated.meta[k] == aia_map_updated2.meta[k]


@pytest.mark.remote_data
def test_fix_pointing_missing_value(aia_171_map):
    # Adjust map to a date we know has missing pointing information
    aia_171_map.meta['date-obs'] = '2010-09-30T05:51:48.344'
    with pytest.warns(AiapyUserWarning, match='Missing value in pointing table'):
        aia_171_map_updated = update_pointing(aia_171_map)
    assert aia_171_map.meta['crpix1'] == aia_171_map_updated.meta['crpix1']
    assert aia_171_map.meta['crpix2'] == aia_171_map_updated.meta['crpix2']
