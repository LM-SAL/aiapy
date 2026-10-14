import copy

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.map.mapbase import PixelPair

from aiapy.calibrate import fetch_spikes, respike
from aiapy.tests.helpers import make_aia_map
from aiapy.utils import AIApyUserWarning


@pytest.fixture
def synthetic_map():
    return make_aia_map(np.zeros((8, 8)), 193, level=1, nspikes=2)


# Module scoped so each JSOC query is made once
@pytest.fixture(scope="module")
def respiked_map(aia_193_level1_map):
    return respike(aia_193_level1_map)


@pytest.fixture(scope="module")
def spikes(aia_193_level1_map):
    return fetch_spikes(aia_193_level1_map)


@pytest.mark.remote_data
def test_respike(respiked_map, spikes) -> None:
    coords, values = spikes
    for x, y, v in zip(coords.x.value, coords.y.value, values, strict=True):
        assert v == respiked_map.data[int(y), int(x)]


def test_respike_supplied_spikes(synthetic_map) -> None:
    # WCS pixel positions can be slightly off integers; they must round to the nearest pixel.
    coords = PixelPair(x=[1, 3.0000001] * u.pix, y=[2, 4.9999999] * u.pix)
    respiked = respike(synthetic_map, spikes=(coords, np.array([10.0, 20.0])))
    expected = np.zeros((8, 8))
    expected[2, 1] = 10.0
    expected[5, 3] = 20.0
    np.testing.assert_array_equal(respiked.data, expected)
    np.testing.assert_array_equal(synthetic_map.data, 0)
    assert respiked.meta["lvl_num"] == 0.5
    assert respiked.meta["nspikes"] == 0


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


@pytest.mark.parametrize(
    ("key", "value", "error", "match"),
    [
        ("lvl_num", 1.5, ValueError, "Can only apply respike procedure to level 1 data"),
        ("nspikes", 0, ValueError, "No spikes were present in the level 0 data."),
        ("instrume", "not AIA", TypeError, "Input must be an AIAMap."),
    ],
)
def test_exceptions(synthetic_map, key, value, error, match) -> None:
    new_meta = copy.deepcopy(synthetic_map.meta)
    new_meta[key] = value
    with pytest.raises(error, match=match):
        respike(sunpy.map.Map(synthetic_map.data, new_meta))


def test_resample_warning(synthetic_map) -> None:
    resampled = synthetic_map.resample((4, 4) * u.pixel)
    spikes = (PixelPair(x=[1] * u.pix, y=[1] * u.pix), np.array([10.0]))
    with pytest.warns(AIApyUserWarning, match="is significantly different from the expected level 1 plate scale"):
        respike(resampled, spikes=spikes)


@pytest.mark.remote_data
def test_fetch_spikes(aia_193_level1_map, spikes) -> None:
    n_spikes = aia_193_level1_map.meta["nspikes"]
    pixels, values = spikes
    assert isinstance(pixels, PixelPair)
    assert all(c.size == n_spikes for c in pixels)
    assert values.size == n_spikes
    coords, values = fetch_spikes(aia_193_level1_map, as_coords=True)
    assert isinstance(coords, SkyCoord)
    assert len(coords) == n_spikes
    assert values.size == n_spikes
