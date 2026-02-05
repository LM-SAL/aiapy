import copy

import numpy as np
import pytest

import astropy.units as u
from astropy.table import QTable
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time, TimeDelta

from aiapy.calibrate import update_pointing
from aiapy.calibrate.utils import get_pointing_table
from aiapy.utils.exceptions import AIApyUserWarning


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
            "A_193_X0",
            "A_193_Y0",
            "A_193_INSTROT",
            "A_193_IMSCALE",
        ),
    )


@pytest.mark.remote_data
def test_fix_pointing(aia_193_level1_map, pointing_table) -> None:
    # Smoke test to make sure expected keys are being updated
    keys = ["crpix1", "crpix2", "cdelt1", "cdelt2", "crota2", "x0_mp", "y0_mp"]
    # NOTE: This modification forces CROTA2 and CDELT{1,2} to be updated. Otherwise,
    # these keys are not modified since the existing metadata and the pointing table
    # for this time are the same.
    new_pointing_table = pointing_table.copy()
    new_pointing_table["A_193_INSTROT"] = 0 * u.degree
    new_pointing_table["A_193_IMSCALE"] = 1.2 * u.arcsec / u.pixel
    aia_map_updated = update_pointing(aia_193_level1_map, pointing_table=new_pointing_table)
    for k in keys:
        assert k in aia_map_updated.meta.modified_items


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
def test_update_pointing_accuracy(aia_193_level1_map, pointing_table, t_delt_factor, expected_entry) -> None:
    t_start = pointing_table[0]["T_START"]
    t_delt = pointing_table[0]["T_STOP"] - t_start  # This is nearly always 3 hours
    aia_193_level1_map.meta["T_OBS"] = (t_start + t_delt * t_delt_factor).isot
    aia_map_updated = update_pointing(aia_193_level1_map, pointing_table=pointing_table)
    assert aia_map_updated.reference_pixel.x == pointing_table[expected_entry]["A_193_X0"]
    assert aia_map_updated.reference_pixel.y == pointing_table[expected_entry]["A_193_Y0"]


@pytest.mark.remote_data
def test_update_pointing_submap(aia_193_level1_map, pointing_table) -> None:
    # Tests that submapping and update pointing are commutative
    # NOTE: Purposefully cropping in pixel space as cropping in world space
    # can lead to differences of 1 pixel due to the differing world-to-pixel
    # conversions because of the change in reference pixel.
    blc = (500, 1000) * u.pixel
    trc = (1000, 1200) * u.pixel
    m_submap_pointing_update = update_pointing(
        aia_193_level1_map.submap(blc, top_right=trc), pointing_table=pointing_table
    )
    m_pointing_update_submap = update_pointing(aia_193_level1_map, pointing_table=pointing_table).submap(
        blc, top_right=trc
    )
    assert_quantity_allclose(
        u.Quantity(m_submap_pointing_update.reference_pixel),
        u.Quantity(m_pointing_update_submap.reference_pixel),
        rtol=1e-10,
    )


@pytest.mark.remote_data
def test_update_pointing_resampled_warning(aia_193_level1_map, pointing_table) -> None:
    m = aia_193_level1_map.resample((512, 512) * u.pixel)
    with pytest.warns(AIApyUserWarning, match="Input map has plate scale"):
        update_pointing(m, pointing_table=pointing_table)


@pytest.mark.remote_data
def test_update_pointing_missing_mp_keys(aia_193_level1_map, pointing_table) -> None:
    new_meta = copy.deepcopy(aia_193_level1_map.meta)
    del new_meta["x0_mp"]
    m = aia_193_level1_map._new_instance(aia_193_level1_map.data, new_meta)
    with pytest.warns(AIApyUserWarning, match="x0_mp and/or y0_mp keywords are missing."):
        m_updated = update_pointing(m, pointing_table=pointing_table)
    assert m.reference_pixel.x == m_updated.reference_pixel.x
    assert m.reference_pixel.y == m_updated.reference_pixel.y


@pytest.mark.remote_data
def test_update_pointing_no_entry_raises_exception(aia_193_level1_map, pointing_table) -> None:
    # This tests that an exception is thrown when entry corresponding to
    # T_START <= T_OBS < T_END cannot be found in the pointing table.
    # We explicitly set the T_OBS key
    aia_193_level1_map.meta["T_OBS"] = (aia_193_level1_map.date - 30 * u.year).isot
    with pytest.raises(IndexError, match="No valid entries for"):
        update_pointing(aia_193_level1_map, pointing_table=pointing_table)


def test_fix_pointing_missing_value(aia_193_level1_map, mock_pointing_table) -> None:
    # Adjust map to a date we know has missing pointing information
    aia_193_level1_map.meta["date-obs"] = "2010-09-30T06:51:48.344"
    aia_193_level1_map.meta["t_obs"] = aia_193_level1_map.meta["date-obs"]
    with pytest.warns(AIApyUserWarning, match="Missing value in pointing table"):
        aia_171_map_updated = update_pointing(aia_193_level1_map, pointing_table=mock_pointing_table)
    assert aia_193_level1_map.meta["crpix1"] == aia_171_map_updated.meta["crpix1"]
    assert aia_193_level1_map.meta["crpix2"] == aia_171_map_updated.meta["crpix2"]
