import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.table import QTable
from astropy.time import Time, TimeDelta

from aiapy.calibrate import update_pointing
from aiapy.calibrate.util import get_pointing_table
from aiapy.util.exceptions import AIApyUserWarning


@pytest.fixture
def pointing_table():
    return get_pointing_table("lmsal")


@pytest.fixture
def mock_pointing_table():
    return QTable(
        [
            Time("2010-09-29 21:00:00") + TimeDelta(3 * u.hour) * np.linspace(1, 8, 8),
            Time("2010-09-30 00:00:00") + TimeDelta(3 * u.hour) * np.linspace(1, 8, 8),
            [np.nan * u.pix] * 8,
            [np.nan * u.pix] * 8,
            [0.019327 * u.degree] * 8,
            [0.019327 * u.arcsec / u.pix] * 8,
        ],
        names=(
            "T_START",
            "T_STOP",
            "A_171_X0",
            "A_171_Y0",
            "A_171_INSTROT",
            "A_171_IMSCALE",
        ),
    )


@pytest.mark.remote_data
def test_fix_pointing(aia_171_map, pointing_table) -> None:
    keys = ["CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "CROTA2"]
    # Remove keys to at least test that they get set
    for k in keys:
        aia_171_map.meta.pop(k)
    aia_map_updated = update_pointing(
        aia_171_map,
        pointing_table=pointing_table,
    )
    # FIXME: how do we check these values are accurate?
    assert all(k in aia_map_updated.meta for k in keys)
    # Check the case where we have specified the pointing
    # table ahead of time
    aia_map_updated2 = update_pointing(aia_171_map, pointing_table=pointing_table)
    for k in keys:
        assert aia_map_updated.meta[k] == aia_map_updated2.meta[k]


@pytest.mark.remote_data
@pytest.mark.parametrize(
    ("t_delt_factor", "expected_entry"),
    [
        # T_OBS = T_START[i] chooses the i-th entry
        (0, 0),
        # T_OBS = T_START[i] + epsilon chooses the i-th entry
        (0.001, 0),
        # T_OBS = (T_START[i] + T_STOP[i])/2 + epsilon chooses the i-th entry
        (0.501, 0),
        # T_OBS = T_STOP[i] chooses the i+1-th entry
        (1, 1),
        # T_OBS = T_STOP[i] + epsilon chooses the i+1-th entry
        (1.001, 1),
    ],
)
def test_update_pointing_accuracy(aia_171_map, pointing_table, t_delt_factor, expected_entry) -> None:
    t_start = pointing_table[0]["T_START"]
    t_delt = pointing_table[0]["T_STOP"] - t_start  # This is nearly always 3 hours
    aia_171_map.meta["T_OBS"] = (t_start + t_delt * t_delt_factor).isot
    aia_map_updated = update_pointing(aia_171_map, pointing_table=pointing_table)
    assert aia_map_updated.reference_pixel.x == pointing_table[expected_entry]["A_171_X0"]
    assert aia_map_updated.reference_pixel.y == pointing_table[expected_entry]["A_171_Y0"]


@pytest.mark.remote_data
def test_update_pointing_submap_raises_exception(aia_171_map, pointing_table) -> None:
    m = aia_171_map.submap(
        SkyCoord(0, 0, unit="arcsec", frame=aia_171_map.coordinate_frame),
        top_right=aia_171_map.top_right_coord,
    )
    with pytest.raises(ValueError, match="Input must be a full disk image."):
        update_pointing(m, pointing_table=pointing_table)


@pytest.mark.remote_data
def test_update_pointing_resampled_raises_exception(aia_171_map, pointing_table) -> None:
    m = aia_171_map.resample((512, 512) * u.pixel)
    with pytest.raises(ValueError, match="Input must be at the full resolution"):
        update_pointing(m, pointing_table=pointing_table)


@pytest.mark.remote_data
def test_update_pointing_no_entry_raises_exception(aia_171_map, pointing_table) -> None:
    # This tests that an exception is thrown when entry corresponding to
    # T_START <= T_OBS < T_END cannot be found in the pointing table.
    # We explicitly set the T_OBS key
    aia_171_map.meta["T_OBS"] = (aia_171_map.date - 1000 * u.day).isot
    with pytest.raises(IndexError, match="No valid entries for"):
        update_pointing(aia_171_map, pointing_table=pointing_table)


def test_fix_pointing_missing_value(aia_171_map, mock_pointing_table) -> None:
    # Adjust map to a date we know has missing pointing information
    aia_171_map.meta["date-obs"] = "2010-09-30T06:51:48.344"
    aia_171_map.meta["t_obs"] = aia_171_map.meta["date-obs"]
    with pytest.warns(AIApyUserWarning, match="Missing value in pointing table"):
        aia_171_map_updated = update_pointing(aia_171_map, pointing_table=mock_pointing_table)
    assert aia_171_map.meta["crpix1"] == aia_171_map_updated.meta["crpix1"]
    assert aia_171_map.meta["crpix2"] == aia_171_map_updated.meta["crpix2"]
