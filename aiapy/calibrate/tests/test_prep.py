import copy

import numpy as np
import pytest

import astropy.time
import astropy.units as u
from astropy.io.fits.verify import VerifyWarning

import sunpy.data.test
from sunpy.map import Map

from aiapy.calibrate import correct_degradation, degradation, register
from aiapy.calibrate.util import get_correction_table
from aiapy.tests.data import get_test_filepath
from aiapy.util import AIApyUserWarning


@pytest.fixture
def lvl_15_map(aia_171_map):
    return register(aia_171_map)


@pytest.fixture
def non_sdo_map():
    return Map(sunpy.data.test.get_test_filepath("hsi_image_20101016_191218.fits"))


def test_register(aia_171_map, lvl_15_map) -> None:
    """
    Test that header info for the map has been correctly updated after the map
    has been scaled to 0.6 arcsec / pixel and aligned with solar north.
    """
    # TODO: Check all of these for Map attributes and .meta values?
    # Check array shape - We cut off two pixels on each side for kicks
    # Due to fixes in sunpy 3.1.6, the shape is different
    # See https://github.com/sunpy/sunpy/pull/5803
    assert lvl_15_map.data.shape == (4094, 4094) != aia_171_map.data.shape
    # Check crpix values
    assert lvl_15_map.meta["crpix1"] == lvl_15_map.data.shape[1] / 2.0 + 0.5
    assert lvl_15_map.meta["crpix2"] == lvl_15_map.data.shape[0] / 2.0 + 0.5
    # Check cdelt values
    assert lvl_15_map.meta["cdelt1"] / 0.6 == int(lvl_15_map.meta["cdelt1"] / 0.6)
    assert lvl_15_map.meta["cdelt2"] / 0.6 == int(lvl_15_map.meta["cdelt2"] / 0.6)
    # Check rotation value, I am assuming that the inaccuracy in
    # the CROTA -> PCi_j matrix is causing the inaccuracy here
    np.testing.assert_allclose(lvl_15_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert lvl_15_map.meta["lvl_num"] == 1.5


def test_register_filesave(lvl_15_map, tmp_path) -> None:
    """
    Test that adjusted header values are still correct after saving the map and
    reloading it.
    """
    filename = tmp_path / "test_register_filesave.fits"
    with pytest.warns(VerifyWarning, match="The 'BLANK' keyword is only applicable to integer data"):
        lvl_15_map.save(str(filename), overwrite=True)
    load_map = Map(str(filename))
    # Check crpix values
    assert load_map.meta["crpix1"] == lvl_15_map.data.shape[1] / 2.0 + 0.5
    assert load_map.meta["crpix2"] == lvl_15_map.data.shape[0] / 2.0 + 0.5
    # Check cdelt values
    assert load_map.meta["cdelt1"] / 0.6 == int(load_map.meta["cdelt1"] / 0.6)
    assert load_map.meta["cdelt2"] / 0.6 == int(load_map.meta["cdelt2"] / 0.6)
    # Check rotation value
    np.testing.assert_allclose(lvl_15_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert load_map.meta["lvl_num"] == 1.5


def test_register_unsupported_maps(aia_171_map, non_sdo_map) -> None:
    """
    Make sure we raise an error when an unsupported map is passed in.
    """
    # A submap
    original_cutout = aia_171_map.submap(aia_171_map.center, top_right=aia_171_map.top_right_coord)
    with pytest.raises(ValueError, match="Input must be a full disk image."):
        register(original_cutout)
    # A Map besides AIA or HMI
    with pytest.raises(TypeError, match="Input must be an AIAMap"):
        register(non_sdo_map)


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
        pytest.param("jsoc", marks=pytest.mark.remote_data),
        pytest.param("SsW", marks=pytest.mark.remote_data),
        get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
        str(get_test_filepath("aia_V8_20171210_050627_response_table.txt")),
    ],
)
def test_correct_degradation(aia_171_map, source) -> None:
    correction_table = get_correction_table(source=source)
    original_corrected = correct_degradation(
        aia_171_map,
        correction_table=correction_table,
    )
    d = degradation(
        aia_171_map.wavelength,
        aia_171_map.date,
        correction_table=correction_table,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        uncorrected_over_corrected = aia_171_map.data / original_corrected.data
    # If intensity is zero, ratio will be NaN/infinite
    i_valid = aia_171_map.data > 0.0
    assert np.allclose(uncorrected_over_corrected[i_valid], d)


@pytest.mark.parametrize(
    ("source", "time_correction_truth"),
    [
        pytest.param(
            "SSW",
            0.9031773242843387 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            "JSOC",
            0.86288462 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        (
            get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
            0.7667108920899671 * u.dimensionless_unscaled,
        ),
    ],
)
def test_degradation(source, time_correction_truth) -> None:
    # NOTE: this just tests an expected result from aiapy, not necessarily an
    # absolutely correct result. It was calculated for the above time and
    # the specific correction table file.
    # TODO: Test this over multiple wavelengths
    correction_table = get_correction_table(source=source)
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    time_correction = degradation(
        94 * u.angstrom,
        obstime,
        correction_table=correction_table,
    )
    assert u.allclose(time_correction, time_correction_truth, atol=1e-3)


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
    assert u.allclose(time_correction, result, atol=1e-3)


@pytest.mark.remote_data
def test_degradation_4500_missing() -> None:
    # 4500 has a max version of 3, so by default it will error
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    with pytest.raises(
        ValueError,
        match="Correction table does not contain calibration for 4500.0 Angstrom. Max version is 3.",
    ):
        degradation(4500 * u.angstrom, obstime, correction_table=get_correction_table("SSW"))


@pytest.mark.remote_data
def test_degradation_4500_jsoc() -> None:
    # 4500 has a max version of 3, so by default it will error
    # and it is missing from the SSW files but not the JSOC
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    correction = degradation(4500 * u.angstrom, obstime, correction_table=get_correction_table("jsoc"))
    assert u.allclose(correction, 1.0 * u.dimensionless_unscaled)


def test_degradation_time_array() -> None:
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    obstime = obstime + np.linspace(0, 1, 100) * u.year
    correction_table = get_correction_table(get_test_filepath("aia_V8_20171210_050627_response_table.txt"))
    time_correction = degradation(
        94 * u.angstrom,
        obstime,
        correction_table=correction_table,
    )
    assert time_correction.shape == obstime.shape
    for o, tc in zip(obstime, time_correction, strict=True):
        assert tc == degradation(94 * u.angstrom, o, correction_table=correction_table)


def test_register_cupy(aia_171_map) -> None:
    pytest.importorskip("cupy")
    cupy_map = register(aia_171_map, method="cupy")
    scipy_map = register(aia_171_map, method="scipy")
    assert np.allclose(cupy_map.data, scipy_map.data)
