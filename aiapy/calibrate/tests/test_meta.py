import numpy as np
import pytest

import astropy.units as u
from astropy.table import QTable
from astropy.tests.helper import assert_quantity_allclose

from aiapy.calibrate import update_pointing
from aiapy.tests.helpers import make_aia_map
from aiapy.utils.exceptions import AIApyUserWarning


@pytest.fixture
def synthetic_map():
    # The reference pixel (CRPIX - 1) matches the Sun center recorded in X0_MP/Y0_MP.
    return make_aia_map(np.zeros((8, 8)), 193, level=1, x0_mp=3.5, y0_mp=3.5, sat_rot=0.1)


@pytest.fixture
def synthetic_pointing_table(synthetic_map):
    # Two consecutive 3 hour entries; the map's T_OBS falls in the first one.
    t_start = synthetic_map.reference_date + [-1, 2] * u.h
    return QTable(
        [
            t_start,
            t_start + 3 * u.h,
            [4.5, 5.5] * u.pix,
            [2.5, 1.5] * u.pix,
            [0.02, 0.03] * u.deg,
            [0.6, 0.61] * u.arcsec / u.pix,
        ],
        names=("T_START", "T_STOP", "A_193_X0", "A_193_Y0", "A_193_INSTROT", "A_193_IMSCALE"),
    )


def test_update_pointing(synthetic_map, synthetic_pointing_table) -> None:
    updated = update_pointing(synthetic_map, pointing_table=synthetic_pointing_table)
    # CRPIX moves with the Sun center; CROTA2 is INSTROT + SAT_ROT.
    assert updated.meta["crpix1"] == synthetic_map.meta["crpix1"] + 1
    assert updated.meta["crpix2"] == synthetic_map.meta["crpix2"] - 1
    assert updated.meta["x0_mp"] == 4.5
    assert updated.meta["y0_mp"] == 2.5
    assert_quantity_allclose(updated.scale, [0.6, 0.6] * u.arcsec / u.pix)
    assert updated.meta["crota2"] == pytest.approx(0.12)
    np.testing.assert_array_equal(updated.data, synthetic_map.data)


@pytest.mark.parametrize(
    ("t_delta_factor", "expected_entry"),
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
def test_update_pointing_accuracy(synthetic_map, synthetic_pointing_table, t_delta_factor, expected_entry) -> None:
    t_start = synthetic_pointing_table[0]["T_START"]
    t_delta = synthetic_pointing_table[0]["T_STOP"] - t_start
    synthetic_map.meta["t_obs"] = (t_start + t_delta * t_delta_factor).isot
    updated = update_pointing(synthetic_map, pointing_table=synthetic_pointing_table)
    assert updated.reference_pixel.x == synthetic_pointing_table[expected_entry]["A_193_X0"]
    assert updated.reference_pixel.y == synthetic_pointing_table[expected_entry]["A_193_Y0"]


def test_update_pointing_missing_value(synthetic_map, synthetic_pointing_table) -> None:
    synthetic_pointing_table["A_193_X0"][0] = np.nan * u.pix
    with pytest.warns(AIApyUserWarning, match="Missing value in pointing table"):
        updated = update_pointing(synthetic_map, pointing_table=synthetic_pointing_table)
    # Only the keys derived from the missing value are left unchanged.
    assert updated.meta["crpix1"] == synthetic_map.meta["crpix1"]
    assert updated.meta["x0_mp"] == synthetic_map.meta["x0_mp"]
    assert updated.meta["crpix2"] == synthetic_map.meta["crpix2"] - 1


def test_update_pointing_missing_mp_keys(synthetic_map, synthetic_pointing_table) -> None:
    del synthetic_map.meta["x0_mp"]
    with pytest.warns(AIApyUserWarning, match="x0_mp and/or y0_mp keywords are missing."):
        updated = update_pointing(synthetic_map, pointing_table=synthetic_pointing_table)
    assert updated.reference_pixel.x == synthetic_map.reference_pixel.x
    assert updated.reference_pixel.y == synthetic_map.reference_pixel.y


def test_update_pointing_no_entry_raises_exception(synthetic_map, synthetic_pointing_table) -> None:
    synthetic_map.meta["t_obs"] = (synthetic_map.reference_date - 30 * u.year).isot
    with pytest.raises(IndexError, match="No valid entries for"):
        update_pointing(synthetic_map, pointing_table=synthetic_pointing_table)


def test_update_pointing_resampled_warning(synthetic_map, synthetic_pointing_table) -> None:
    resampled = synthetic_map.resample((4, 4) * u.pixel)
    with pytest.warns(AIApyUserWarning, match="Input map has plate scale"):
        update_pointing(resampled, pointing_table=synthetic_pointing_table)


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
