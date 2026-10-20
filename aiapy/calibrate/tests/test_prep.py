import copy
from unittest.mock import Mock

import numpy as np
import pytest

import astropy.time
import astropy.units as u
from astropy.io.fits.verify import VerifyWarning
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs.utils import pixel_to_pixel

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.sdo import HMIMap

from aiapy.calibrate import correct_degradation, degradation, register
from aiapy.calibrate.utils import get_correction_table
from aiapy.tests.data import get_test_filepath
from aiapy.tests.helpers import make_aia_map
from aiapy.utils import AIApyUserWarning, detector_dimensions


@pytest.fixture(scope="module")
def lvl_15_map():
    smap = Map(sunpy.data.test.get_test_filepath("aia_171_level1.fits"))
    return register(smap.resample(detector_dimensions()))


def time_keywords(meta):
    return {k: v for k, v in meta.items() if k.startswith(("date", "mjd", "t_obs", "timesys"))}


def assert_contains_input(smap, registered):
    """
    Every input pixel corner lands inside the registered output.
    """
    ny, nx = smap.data.shape
    x, y = pixel_to_pixel(
        smap.wcs,
        registered.wcs,
        np.array([-0.5, nx - 0.5, nx - 0.5, -0.5]),
        np.array([-0.5, -0.5, ny - 0.5, ny - 0.5]),
    )
    assert np.all(x >= -0.5 - 1e-7)
    assert np.all(x <= registered.data.shape[1] - 0.5 + 1e-7)
    assert np.all(y >= -0.5 - 1e-7)
    assert np.all(y <= registered.data.shape[0] - 0.5 + 1e-7)


def test_register(lvl_15_map, tmp_path) -> None:
    """
    Check the full-disk level 1.5 grid and calibration metadata, before and after saving.
    """
    filename = tmp_path / "registered.fits"
    with pytest.warns(VerifyWarning, match="The 'BLANK' keyword is only applicable to integer data"):
        lvl_15_map.save(filename)
    for smap in (lvl_15_map, Map(filename)):
        assert smap.data.shape == tuple(detector_dimensions().to_value(u.pix))
        assert_quantity_allclose(smap.reference_pixel, (u.Quantity(smap.dimensions) - 1 * u.pix) / 2)
        assert_quantity_allclose(smap.scale, [0.6, 0.6] * u.arcsec / u.pix)
        assert smap.meta["cunit1"] == smap.meta["cunit2"] == "arcsec"
        np.testing.assert_allclose(smap.rotation_matrix, np.identity(2))
        assert smap.meta["lvl_num"] == 1.5
        assert smap.meta["bitpix"] == -64


def test_register_aia_channels_common_grid(small_aia_map) -> None:
    """
    Different native AIA scales register onto the same pixel grid.
    """
    registered = []
    for wavelength, scale in ((94, 0.599), (171, 0.6), (193, 0.601)):
        smap = make_aia_map(
            small_aia_map.data,
            wavelength,
            scale=scale,
            rotation=23,
            coordinate=small_aia_map.reference_coordinate,
            unit=small_aia_map.unit,
            level=1,
        )
        result = register(smap)
        assert_quantity_allclose(result.scale, [0.6, 0.6] * u.arcsec / u.pix)
        registered.append(result)
    for result in registered[1:]:
        x, y = pixel_to_pixel(result.wcs, registered[0].wcs, np.array([0, 1]), np.array([0, 1]))
        np.testing.assert_allclose([x, y], np.rint([x, y]), atol=1e-7)
        np.testing.assert_allclose([np.diff(x), np.diff(y)], 1, atol=1e-7)


def test_register_hmi_on_aia_grid(hmi_map) -> None:
    """
    HMI registers onto the AIA 0.6 arcsec grid so HMI and AIA images overlay pixel for pixel.
    """
    with pytest.warns(AIApyUserWarning, match="Image registration should only be applied to level 1 data"):
        registered, footprint = register(hmi_map, parallel=False, return_footprint=True)
    assert isinstance(registered, HMIMap)
    assert_quantity_allclose(registered.scale, [0.6, 0.6] * u.arcsec / u.pix)
    # Disk center falls on the same pixel grid as a registered full-disk AIA image.
    offset = u.Quantity(registered.reference_pixel).to_value(u.pix) - (detector_dimensions().to_value(u.pix) - 1) / 2
    np.testing.assert_allclose(offset, np.rint(offset), atol=1e-7)
    assert registered.unit == hmi_map.unit
    assert registered.date == hmi_map.date
    assert registered.reference_date == hmi_map.reference_date
    np.testing.assert_allclose(registered.rotation_matrix, np.identity(2))
    assert np.count_nonzero(footprint) > 0
    assert_contains_input(hmi_map, registered)


@pytest.mark.parametrize("matrix", ["pc", "cd", "crota"])
def test_register_wcs_metadata(small_aia_map, matrix) -> None:
    meta = small_aia_map.meta.copy()
    rotation = small_aia_map.rotation_matrix
    for key in ("pc1_1", "pc1_2", "pc2_1", "pc2_2"):
        meta.pop(key, None)
    if matrix == "crota":
        meta["crota2"] = 23
    else:
        for i in range(2):
            for j in range(2):
                meta[f"{matrix}{i + 1}_{j + 1}"] = rotation[i, j] * (0.6 if matrix == "cd" else 1)
    smap = Map(small_aia_map.data, meta)
    original_meta = smap.meta.copy()
    original_data = smap.data.copy()
    registered = register(smap)
    np.testing.assert_allclose(registered.rotation_matrix, np.identity(2))
    assert_quantity_allclose(registered.scale, [0.6, 0.6] * u.arcsec / u.pix)
    assert registered.meta["naxis1"] == registered.data.shape[1]
    assert registered.meta["naxis2"] == registered.data.shape[0]
    assert smap.meta == original_meta
    np.testing.assert_array_equal(smap.data, original_data)


def test_register_resampled_full_disk(aia_171_test_map) -> None:
    smap = aia_171_test_map
    registered = register(smap)
    assert registered.data.shape == smap.data.shape == (128, 128)
    assert_quantity_allclose(registered.scale, [19.2, 19.2] * u.arcsec / u.pix)
    assert_quantity_allclose(registered.reference_pixel, [63.5, 63.5] * u.pix)
    assert abs(registered.bottom_left_coord.Tx) > 1200 * u.arcsec


@pytest.mark.parametrize("date_keyword", ["date-obs", "date_obs"])
@pytest.mark.parametrize("timesys", ["UTC", "TAI"])
@pytest.mark.parametrize("reference_keyword", ["date-avg", "t_obs"])
def test_register_observation_time(small_aia_map, tmp_path, date_keyword, timesys, reference_keyword) -> None:
    meta = small_aia_map.meta.copy()
    meta.pop("date-obs")
    # The fixture carries reference-time keywords; drop the candidates so this test
    # controls exactly which one (date-avg or t_obs) drives the reference date.
    meta.pop("date-avg", None)
    meta.pop("t_obs", None)
    meta[date_keyword] = "2020-01-01T00:00:00.123456"
    meta["timesys"] = timesys
    reference_date = astropy.time.Time("2020-01-01T00:00:01", scale=timesys.lower())
    meta[reference_keyword] = (
        "2020.01.01_00:00:01_TAI" if reference_keyword == "t_obs" and timesys == "TAI" else reference_date.isot
    )
    meta["mjd-avg"] = reference_date.mjd
    smap = Map(small_aia_map.data, meta)
    assert smap.reference_date == reference_date
    assert smap.date == astropy.time.Time("2020-01-01T00:00:00.123456", scale=timesys.lower())
    registered = register(smap)
    assert time_keywords(registered.meta) == time_keywords(smap.meta)
    filename = tmp_path / "registered.fits"
    registered.save(filename)
    for result in (registered, Map(filename)):
        assert result.date == smap.date
        assert result.reference_date == smap.reference_date


def test_register_solar_radius(aia_171_test_map) -> None:
    smap = aia_171_test_map
    smap.meta["rsun_obs"] = 960
    registered = register(smap)
    assert_quantity_allclose(registered.rsun_obs, smap.rsun_obs)
    assert registered.meta["r_sun"] == pytest.approx(50)
    assert registered.meta["r_sun"] == pytest.approx((registered.rsun_obs / registered.scale.axis1).to_value(u.pix))


def test_register_real_input_metadata(aia_171_test_map) -> None:
    """
    A real AIA level-1 header (T_OBS with a trailing Z, no TIMESYS or MJD keywords) keeps its times.
    """
    registered = register(aia_171_test_map)
    assert time_keywords(registered.meta) == time_keywords(aia_171_test_map.meta)


@pytest.mark.parametrize("algorithm", ["interpolation", "adaptive", "exact"])
@pytest.mark.parametrize("missing", [None, -100])
def test_register_footprint_and_missing(small_aia_map, algorithm, missing) -> None:
    registered, footprint = register(small_aia_map, algorithm=algorithm, missing=missing, return_footprint=True)
    assert footprint.shape == registered.data.shape
    assert np.any(footprint > 0)
    assert np.any(footprint == 0)
    assert np.all(np.isfinite(registered.data))
    np.testing.assert_array_equal(registered.data[footprint == 0], small_aia_map.min() if missing is None else missing)


def test_register_no_overlap(small_aia_map) -> None:
    # Keep a target extent far from the cutout instead of automatically expanding it.
    registered, footprint = register(small_aia_map, auto_extent=None, shape_out=(32, 40), return_footprint=True)
    np.testing.assert_array_equal(footprint, 0)
    np.testing.assert_array_equal(registered.data, small_aia_map.min())


@pytest.mark.parametrize("order", [1, 3])
@pytest.mark.parametrize("missing", [None, -100, np.nan])
def test_register_nan_input(small_aia_map, order, missing) -> None:
    data = small_aia_map.data.astype(float)
    data[15, 20] = np.nan
    original_data = data.copy()
    smap = small_aia_map._new_instance(data, small_aia_map.meta)
    registered, footprint = register(smap, order=order, missing=missing, return_footprint=True)
    assert np.count_nonzero(footprint) > data.size / 2
    np.testing.assert_array_equal(registered.data[footprint == 0], smap.min() if missing is None else missing)
    assert np.all(np.isfinite(registered.data[footprint > 0]))
    assert np.max(registered.data[footprint > 0]) <= smap.max()
    assert np.min(registered.data[footprint > 0]) >= smap.min()
    np.testing.assert_array_equal(smap.data, original_data)


@pytest.mark.parametrize(("available_cpus", "expected_parallel"), [(None, 1), (8, 8)])
def test_register_parallel_defaults(small_aia_map, monkeypatch, available_cpus, expected_parallel) -> None:
    monkeypatch.setattr("aiapy.calibrate.prep.os.process_cpu_count", lambda: available_cpus)
    large_map = small_aia_map._new_instance(np.broadcast_to(1, (4096, 4096)), small_aia_map.meta)
    reprojection = Mock(return_value=small_aia_map)
    monkeypatch.setattr(large_map, "reproject_to", reprojection)
    register(large_map)
    assert reprojection.call_args.kwargs["auto_extent"] is None
    assert reprojection.call_args.kwargs["block_size"] == (1024, 1024)
    assert reprojection.call_args.kwargs["parallel"] == expected_parallel
    assert reprojection.call_args.kwargs["order"] == 3
    assert reprojection.call_args.kwargs["roundtrip_coords"] is False
    register(large_map, parallel=False, block_size=(512, 512), order=1)
    assert reprojection.call_args.kwargs["block_size"] == (512, 512)
    assert reprojection.call_args.kwargs["parallel"] is False
    assert reprojection.call_args.kwargs["order"] == 1


def test_register_submap(lvl_15_map, aia_171_submap) -> None:
    """
    A registered cutout keeps its field of view and matches the registered full-disk grid and data.
    """
    submap_lvl_15 = register(aia_171_submap)
    assert_contains_input(aia_171_submap, submap_lvl_15)
    offset = np.array(pixel_to_pixel(submap_lvl_15.wcs, lvl_15_map.wcs, 0, 0))
    np.testing.assert_allclose(offset, np.rint(offset), atol=1e-7)
    x0, y0 = np.rint(offset).astype(int)
    ny, nx = submap_lvl_15.data.shape
    rows, columns = slice(ny // 2 - 16, ny // 2 + 16), slice(nx // 2 - 16, nx // 2 + 16)
    full_rows = slice(y0 + rows.start, y0 + rows.stop)
    full_columns = slice(x0 + columns.start, x0 + columns.stop)
    np.testing.assert_allclose(submap_lvl_15.data[rows, columns], lvl_15_map.data[full_rows, full_columns], rtol=1e-7)


def test_register_level_15(small_aia_map) -> None:
    meta = copy.deepcopy(small_aia_map.meta)
    meta["lvl_num"] = 1.5
    with pytest.warns(AIApyUserWarning, match="Image registration should only be applied to level 1 data"):
        register(small_aia_map._new_instance(small_aia_map.data, meta))
    # A missing processing level also warns.
    del meta["lvl_num"]
    with pytest.warns(AIApyUserWarning, match="Image registration should only be applied to level 1 data"):
        register(small_aia_map._new_instance(small_aia_map.data, meta))


@pytest.mark.parametrize(
    ("source"),
    [
        pytest.param(None, marks=pytest.mark.remote_data),
        get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
        str(get_test_filepath("aia_V8_20171210_050627_response_table.txt")),
    ],
)
def test_correct_degradation(aia_171_test_map, source) -> None:
    correction_table = get_correction_table() if source is None else get_correction_table(source=source)
    calibration_version = np.max(correction_table["VER_NUM"])
    original_corrected = correct_degradation(
        aia_171_test_map,
        correction_table=correction_table,
        calibration_version=calibration_version,
    )
    d = degradation(
        aia_171_test_map.wavelength,
        aia_171_test_map.date,
        correction_table=correction_table,
        calibration_version=calibration_version,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        uncorrected_over_corrected = aia_171_test_map.data / original_corrected.data
    # If intensity is zero, ratio will be NaN/infinite
    i_valid = aia_171_test_map.data > 0.0
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


@pytest.mark.remote_data
@pytest.mark.parametrize(
    ("wavelength", "result"),
    [
        (1600, 0.56048306),
        (1700, 0.87895035),
        (171, 0.78689527),
        (304, 0.14925002),
        (211, 0.82998772),
        (193, 0.86181854),
        (335, 0.32689207),
        (94, 0.90317732),
    ],
)
def test_degradation_all_wavelengths(wavelength, result) -> None:
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    time_correction = degradation(
        wavelength * u.angstrom,
        obstime,
        correction_table=get_correction_table("SSW"),
    )
    assert_quantity_allclose(time_correction, result * u.dimensionless_unscaled, atol=1e-3)


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
def test_degradation_4500_jsoc(jsoc_correction_table) -> None:
    # 4500 has a max version of 3, so by default it will error
    # and it is missing from the SSW files but not the JSOC
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    correction = degradation(4500 * u.angstrom, obstime, correction_table=jsoc_correction_table)
    assert_quantity_allclose(correction, 1.0 * u.dimensionless_unscaled)


def test_degradation_time_array() -> None:
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    obstime = obstime + np.linspace(0, 1, 100) * u.year
    correction_table = get_correction_table(get_test_filepath("aia_V8_20171210_050627_response_table.txt"))
    time_correction = degradation(94 * u.angstrom, obstime, correction_table=correction_table, calibration_version=8)
    assert time_correction.shape == obstime.shape
    for o, tc in zip(obstime, time_correction, strict=True):
        assert tc == degradation(94 * u.angstrom, o, correction_table=correction_table, calibration_version=8)
