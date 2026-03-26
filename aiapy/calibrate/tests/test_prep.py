import copy

import numpy as np
import pytest

import astropy.time
import astropy.units as u
from astropy.io.fits.verify import VerifyWarning
from astropy.tests.helper import assert_quantity_allclose

from sunpy.map import Map

from aiapy.calibrate import correct_degradation, degradation, register
from aiapy.calibrate.utils import get_correction_table
from aiapy.tests.data import get_test_filepath
from aiapy.utils import AIApyUserWarning, detector_dimensions


@pytest.fixture
def lvl_15_map(aia_171_map):
    return register(aia_171_map)


def test_register(lvl_15_map) -> None:
    """
    Test that header info for the map has been correctly updated after the map
    has been scaled to 0.6 arcsec / pixel and aligned with solar north.
    """
    np.testing.assert_array_equal(lvl_15_map.data.shape, detector_dimensions().value)
    assert_quantity_allclose(
        u.Quantity(lvl_15_map.reference_pixel), (u.Quantity(lvl_15_map.dimensions) - 1 * u.pix) / 2
    )
    assert_quantity_allclose(lvl_15_map.scale, [0.6, 0.6] * u.arcsec / u.pixel)
    np.testing.assert_allclose(lvl_15_map.rotation_matrix, np.identity(2))
    assert lvl_15_map.meta["lvl_num"] == 1.5
    assert lvl_15_map.meta["bitpix"] == -64


def test_register_submap(aia_171_submap) -> None:
    """
    Test that submap and register are commutative
    """
    submap_lvl_15 = register(aia_171_submap)
    np.testing.assert_array_equal(submap_lvl_15.data.shape, aia_171_submap.data.shape)
    assert_quantity_allclose(u.Quantity(submap_lvl_15.reference_pixel), u.Quantity(aia_171_submap.reference_pixel))
    assert_quantity_allclose(submap_lvl_15.scale, [0.6, 0.6] * u.arcsec / u.pixel)
    np.testing.assert_allclose(submap_lvl_15.rotation_matrix, np.identity(2))
    assert submap_lvl_15.meta["lvl_num"] == 1.5
    assert submap_lvl_15.meta["bitpix"] == -64


def test_register_submap_commutative(lvl_15_map, aia_171_submap) -> None:
    """
    Test that submap and register are commutative
    """
    submap_lvl_15 = register(aia_171_submap)
    lvl_15_submap = lvl_15_map.submap(aia_171_submap.bottom_left_coord, top_right=aia_171_submap.top_right_coord)
    assert_quantity_allclose(u.Quantity(submap_lvl_15.reference_pixel), u.Quantity(lvl_15_submap.reference_pixel))
    assert submap_lvl_15.data.shape == lvl_15_submap.data.shape


def test_register_filesave(lvl_15_map, tmp_path) -> None:
    """
    Test that adjusted header values are still correct after saving the map and
    reloading it.
    """
    filename = tmp_path / "test_register_filesave.fits"
    with pytest.warns(VerifyWarning, match="The 'BLANK' keyword is only applicable to integer data"):
        lvl_15_map.save(str(filename), overwrite=True)
    load_map = Map(str(filename))
    np.testing.assert_array_equal(load_map.data.shape, detector_dimensions().value)
    assert_quantity_allclose(u.Quantity(load_map.reference_pixel), (u.Quantity(load_map.dimensions) - 1 * u.pix) / 2)
    assert_quantity_allclose(load_map.scale, [0.6, 0.6] * u.arcsec / u.pixel)
    np.testing.assert_allclose(load_map.rotation_matrix, np.identity(2))
    np.testing.assert_allclose(load_map.rotation_matrix, np.identity(2))
    assert load_map.meta["lvl_num"] == 1.5
    assert load_map.meta["bitpix"] == -64


def test_register_level_15(lvl_15_map) -> None:
    with pytest.warns(
        AIApyUserWarning,
        match="Image registration should only be applied to level 1 data",
    ):
        register(lvl_15_map)
    new_meta = copy.deepcopy(lvl_15_map.meta)
    # Test case where processing_level is missing and returns None
    del new_meta["lvl_num"]
    with pytest.warns(
        AIApyUserWarning,
        match="Image registration should only be applied to level 1 data",
    ):
        register(lvl_15_map._new_instance(lvl_15_map.data, new_meta))


@pytest.mark.parametrize(
    ("source"),
    [
        pytest.param(None, marks=pytest.mark.remote_data),
        # We test different casings to make sure it's case insensitive
        pytest.param("jsOc", marks=pytest.mark.remote_data),
        pytest.param("SsW", marks=pytest.mark.remote_data),
        get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
        str(get_test_filepath("aia_V8_20171210_050627_response_table.txt")),
    ],
)
def test_correct_degradation(aia_171_map, source) -> None:
    correction_table = get_correction_table() if source is None else get_correction_table(source=source)
    calibration_version = np.max(correction_table["VER_NUM"])
    original_corrected = correct_degradation(
        aia_171_map,
        correction_table=correction_table,
        calibration_version=calibration_version,
    )
    d = degradation(
        aia_171_map.wavelength,
        aia_171_map.date,
        correction_table=correction_table,
        calibration_version=calibration_version,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        uncorrected_over_corrected = aia_171_map.data / original_corrected.data
    # If intensity is zero, ratio will be NaN/infinite
    i_valid = aia_171_map.data > 0.0
    np.testing.assert_allclose(np.mean(uncorrected_over_corrected[i_valid]), d.to_value())


@pytest.mark.parametrize(
    ("source", "time_correction_truth", "calibration_version"),
    [
        pytest.param(
            "SSW",
            0.9031773242843387 * u.dimensionless_unscaled,
            # This is none to ensure that the default value works
            None,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            "JSOC",
            # This should match the value above
            0.9031773242843387 * u.dimensionless_unscaled,
            None,
            marks=pytest.mark.remote_data,
        ),
        (
            get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
            0.7667108920899671 * u.dimensionless_unscaled,
            8,
        ),
    ],
)
def test_degradation(source, time_correction_truth, calibration_version) -> None:
    # NOTE: this just tests an expected result from aiapy, not necessarily an
    # absolutely correct result. It was calculated for the above time and
    # the specific correction table file.
    correction_table = get_correction_table(source=source)
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    kwargs = {}
    if calibration_version is not None:
        kwargs["calibration_version"] = calibration_version
    time_correction = degradation(
        94 * u.angstrom,
        obstime,
        correction_table=correction_table,
        **kwargs,
    )
    assert_quantity_allclose(time_correction, time_correction_truth, atol=1e-3)


@pytest.mark.parametrize(
    ("wavelength", "result"),
    [
        pytest.param(
            1600,
            0.56048306 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            1700,
            0.87895035 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            171,
            0.78689527 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            304,
            0.14925002 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            211,
            0.82998772 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            193,
            0.86181854 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            335,
            0.32689207 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            94,
            0.90317732 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
    ],
)
def test_degradation_all_wavelengths(wavelength, result) -> None:
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    time_correction = degradation(
        wavelength * u.angstrom,
        obstime,
        correction_table=get_correction_table("SSW"),
    )
    assert_quantity_allclose(time_correction, result, atol=1e-3)


@pytest.mark.remote_data
def test_degradation_4500_missing() -> None:
    # 4500 has a max version of 3, so by default it will error
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    with pytest.raises(
        ValueError,
        match=r"Correction table does not contain calibration for 4500 Angstrom and version 3.",
    ):
        degradation(4500 * u.angstrom, obstime, correction_table=get_correction_table("SSW"))


@pytest.mark.remote_data
def test_degradation_4500_jsoc() -> None:
    # 4500 has a max version of 3, so by default it will error
    # and it is missing from the SSW files but not the JSOC
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    correction = degradation(4500 * u.angstrom, obstime, correction_table=get_correction_table("jsoc"))
    assert_quantity_allclose(correction, 1.0 * u.dimensionless_unscaled)


def test_degradation_time_array() -> None:
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    obstime = obstime + np.linspace(0, 1, 100) * u.year
    correction_table = get_correction_table(get_test_filepath("aia_V8_20171210_050627_response_table.txt"))
    time_correction = degradation(94 * u.angstrom, obstime, correction_table=correction_table, calibration_version=8)
    assert time_correction.shape == obstime.shape
    for o, tc in zip(obstime, time_correction, strict=True):
        assert tc == degradation(94 * u.angstrom, o, correction_table=correction_table, calibration_version=8)
