import copy

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.map.mapbase import PixelPair

from aiapy.calibrate import fetch_spikes, respike
from aiapy.utils import AIApyUserWarning


@pytest.fixture
def respiked_map(aia_193_level1_map):
    return respike(aia_193_level1_map)


@pytest.fixture
def spikes(aia_193_level1_map):
    return fetch_spikes(aia_193_level1_map)


@pytest.mark.remote_data
def test_respike(respiked_map, spikes) -> None:
    coords, values = spikes
    for x, y, v in zip(coords.x.value, coords.y.value, values, strict=True):
        assert v == respiked_map.data[int(y), int(x)]


@pytest.mark.remote_data
def test_respike_meta(respiked_map) -> None:
    assert respiked_map.meta["lvl_num"] == 0.5
    assert respiked_map.meta["nspikes"] == 0


@pytest.mark.remote_data
def test_fetch_with_prefetched_spikes(aia_193_level1_map, respiked_map, spikes) -> None:
    respiked_map_prefetched = respike(aia_193_level1_map, spikes=spikes)
    np.testing.assert_allclose(respiked_map.data, respiked_map_prefetched.data)


@pytest.mark.remote_data
def test_cutout(respiked_map, aia_193_level1_map) -> None:
    blc = (-500, -500) * u.arcsec
    trc = (500, 500) * u.arcsec
    respiked_map_cutout = respiked_map.submap(
        SkyCoord(*blc, frame=respiked_map.coordinate_frame),
        top_right=SkyCoord(*trc, frame=respiked_map.coordinate_frame),
    )
    cutout_map_respiked = respike(
        aia_193_level1_map.submap(
            SkyCoord(*blc, frame=aia_193_level1_map.coordinate_frame),
            top_right=SkyCoord(*trc, frame=aia_193_level1_map.coordinate_frame),
        ),
    )
    np.testing.assert_allclose(respiked_map_cutout.data, cutout_map_respiked.data)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    ("key", "value", "error", "match"),
    [
        ("lvl_num", 1.5, ValueError, "Can only apply respike procedure to level 1 data"),
        ("nspikes", 0, ValueError, "No spikes were present in the level 0 data."),
        ("instrume", "not AIA", TypeError, "Input must be an AIAMap."),
    ],
)
def test_exceptions(aia_193_level1_map, key, value, error, match) -> None:
    new_meta = copy.deepcopy(aia_193_level1_map.meta)
    new_meta[key] = value
    with pytest.raises(error, match=match):
        respike(sunpy.map.Map(aia_193_level1_map.data, new_meta))


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore::ResourceWarning")
def test_resample_warning(aia_193_level1_map) -> None:
    despiked_map_resample = aia_193_level1_map.resample((512, 512) * u.pixel)
    with pytest.warns(AIApyUserWarning, match="is significantly different from the expected level 1 plate scale"):
        respike(despiked_map_resample)


@pytest.mark.remote_data
@pytest.mark.parametrize(("as_coords", "kind"), [(True, SkyCoord), (False, PixelPair)])
def test_fetch_spikes(aia_193_level1_map, as_coords, kind) -> None:
    n_spikes = aia_193_level1_map.meta["nspikes"]
    coords, values = fetch_spikes(aia_193_level1_map, as_coords=as_coords)
    assert isinstance(coords, kind)
    assert values.size == n_spikes
    if as_coords:
        assert n_spikes == len(coords)
    else:
        assert all(c.size == n_spikes for c in coords)
