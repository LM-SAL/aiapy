import numpy as np
import pytest

import astropy.time
import astropy.units as u
from astropy.tests.helper import assert_quantity_allclose

from aiapy.calibrate.utils import get_correction_table
from aiapy.response import Channel
from aiapy.tests.data import get_test_filepath


# Mark all tests which use this fixture as online
@pytest.fixture(params=[pytest.param(None, marks=pytest.mark.remote_data)])
def channel(request):
    return Channel(94 * u.angstrom, instrument_file=request.param)


@pytest.fixture
def channel_properties():
    # Properties that we need, but are not checked by the ABC
    return [
        "wavelength",
        "primary_mirror_reflectance",
        "secondary_mirror_reflectance",
        "focal_plane_filter_transmittance",
        "entrance_filter_transmittance",
    ]


@pytest.fixture
def required_keys():
    return [
        "wave",
        "primary",
        "secondary",
        "fp_filter",
        "ent_filter",
        "geoarea",
        "ccd",
        "platescale",
        "elecperev",
        "elecperdn",
    ]


def test_channel_data_has_keys(channel, required_keys) -> None:
    assert all(k in channel._data for k in required_keys)


def test_channel_wavelength(channel) -> None:
    assert channel.channel == 94 * u.angstrom
    assert channel.name == "94"


def test_telescope_number(channel) -> None:
    assert channel.telescope_number == 4


def test_invalid_channel() -> None:
    with pytest.raises(ValueError, match=r'channel "1.0 Angstrom" not in list of valid channels'):
        Channel(1 * u.angstrom)


def test_channel_properties(channel, channel_properties) -> None:
    # Test that expected properties are present and are quantities.
    # This does not test correctness
    for p in channel_properties:
        assert isinstance(getattr(channel, p), u.Quantity)


def test_nominal_effective_area(channel) -> None:
    # The effective area without crosstalk matches the one stored in the instrument file
    channel.include_crosstalk = False
    assert_quantity_allclose(channel.effective_area(), channel._data["effarea"] * u.cm**2, rtol=1e-6)


@pytest.mark.parametrize(
    ("source", "eve_correction_truth"),
    [
        pytest.param(
            "SSW",
            0.9548415 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        pytest.param(
            "JSOC",
            0.954845 * u.dimensionless_unscaled,
            marks=pytest.mark.remote_data,
        ),
        (
            get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
            1.01403862 * u.dimensionless_unscaled,
        ),
    ],
)
def test_eve_correction(channel, source, eve_correction_truth) -> None:
    # NOTE: this just tests an expected result from aiapy, not necessarily an
    # absolutely correct result. It was calculated for the above time and
    # the correction parameters in JSOC at the time this code was committed/
    # the specific correction table file.
    # NOTE: If the first two test starts failing, it may be because the
    # correction table parameters have been updated in JSOC.
    # NOTE: This third and fourth values are different from first because the
    # first two are the expected result from JSOC and the third and fourth are
    # results from the correction table file in SSW. The result returned by
    # JSOC are not necessarily the same as those in the correction table files
    # in SSW though they should be close.
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    correction_table = get_correction_table(source=source)
    calibration_version = np.max(correction_table["VER_NUM"])
    eve_correction = channel.eve_correction(
        obstime, correction_table=correction_table, calibration_version=calibration_version
    )
    assert_quantity_allclose(eve_correction, eve_correction_truth)
    # The channel's own correction table and version are the defaults
    channel.correction_table = correction_table
    channel.calibration_version = calibration_version
    assert_quantity_allclose(channel.eve_correction(obstime), eve_correction_truth)


def test_eve_correction_without_obstime(channel) -> None:
    # Without an obstime, the EVE correction is still applied, once
    channel.correction_table = get_test_filepath("aia_V8_20171210_050627_response_table.txt")
    channel.calibration_version = 8
    obstime = astropy.time.Time("2020-01-01T00:00:00")
    eve = channel.eve_correction(obstime)
    assert_quantity_allclose(channel.eve_correction(), eve)
    area, area_obstime = channel.effective_area(), channel.effective_area(obstime)
    channel.include_eve_correction = True
    assert_quantity_allclose(channel.degradation(), eve)
    assert_quantity_allclose(channel.effective_area(), area * eve)
    assert_quantity_allclose(channel.effective_area(obstime), area_obstime * eve)


@pytest.mark.remote_data
def test_eff_wvln_constant(jsoc_correction_table) -> None:
    # eve_correction() without an obstime uses the first epoch, which relies on this
    for table in [jsoc_correction_table, get_correction_table()]:
        for group in table.group_by(["WAVE_STR", "VER_NUM"]).groups:
            assert np.all(group["EFF_WVLN"] == group["EFF_WVLN"][0])


@pytest.mark.parametrize(
    ("obstime", "correction_table", "calibration_version", "include_eve_correction", "include_crosstalk"),
    [
        (None, None, None, False, True),
        (None, None, None, False, False),
        ("2020-01-01", get_test_filepath("aia_V8_20171210_050627_response_table.txt"), 8, False, True),
        ("2020-01-01", get_test_filepath("aia_V8_20171210_050627_response_table.txt"), 8, False, False),
        (astropy.time.Time.now(), get_test_filepath("aia_V8_20171210_050627_response_table.txt"), 8, True, True),
    ],
)
def test_effective_area(
    channel, obstime, correction_table, calibration_version, include_eve_correction, include_crosstalk
) -> None:
    # NOTE: this does not test correctness, but just that the method can
    # be run with the various combinations of inputs. The tests below test the
    # correctness of the output as evaluated by their similarity to those
    # results from SSW.
    channel.include_crosstalk = include_crosstalk
    channel.include_eve_correction = include_eve_correction
    channel.calibration_version = calibration_version
    channel.correction_table = correction_table
    assert channel.effective_area(obstime=obstime).shape == channel.wavelength.shape


def test_effective_area_obstime(channel) -> None:
    channel.correction_table = get_test_filepath("aia_V8_20171210_050627_response_table.txt")
    channel.calibration_version = 8
    obstime = astropy.time.Time("2020-01-01T00:00:00")
    # The effective area at obstime is the pristine effective area scaled by
    # the time-dependent degradation
    assert_quantity_allclose(channel.effective_area(obstime), channel.effective_area() * channel.degradation(obstime))


@pytest.mark.parametrize(
    ("include_crosstalk", "expected"),
    # Reference values at 94 and 304 Å from SSW: aia_get_response(/area,/dn,evenorm=0)
    # with /noblend for no crosstalk, see test_idl.py::test_wavelength_response_consistent
    [
        (True, [0.5724313855171204, 1.0914278902873775e-07]),
        (False, [0.5724313855171204, 2.958600957948072e-15]),
    ],
)
def test_wavelength_response_matches_idl(channel, include_crosstalk, expected) -> None:
    channel.include_crosstalk = include_crosstalk
    r = channel.wavelength_response() / channel.pixel_solid_angle
    indices = [np.argmin(np.abs(channel.wavelength - w)) for w in [94, 304] * u.angstrom]
    assert_quantity_allclose(r[indices], expected * u.cm**2 * u.DN / u.ph, rtol=1e-4)


@pytest.mark.remote_data
@pytest.mark.parametrize(
    ("channel_wavelength", "expected"),
    # Reference values at the channel wavelength and 5000 Å from SSW: aia_get_response(/area,/dn,/uv,evenorm=0)
    [
        (1600 * u.angstrom, [0.002308473689481616, 3.016716737525016e-11]),
        (1700 * u.angstrom, [0.004395448137074709, 1.8126247169902854e-09]),
        (4500 * u.angstrom, [1.6016370864235796e-05, 8.497925563233366e-08]),
    ],
)
def test_fuv_wavelength_response_matches_idl(channel_wavelength, expected) -> None:
    channel = Channel(channel_wavelength)
    r = channel.wavelength_response() / channel.pixel_solid_angle
    indices = [np.argmin(np.abs(channel.wavelength - w)) for w in [channel_wavelength, 5000 * u.angstrom]]
    assert_quantity_allclose(r[indices], expected * u.cm**2 * u.DN / u.ph, rtol=1e-4)


@pytest.mark.remote_data
@pytest.mark.parametrize("channel_wavelength", [1600 * u.angstrom, 1700 * u.angstrom, 4500 * u.angstrom])
def test_fuv_channel(channel_wavelength, channel_properties, required_keys) -> None:
    # There are a few corner cases for the 1600, 1700, and 4500 channels
    channel = Channel(channel_wavelength)
    assert all(k in channel._data for k in required_keys)
    for p in channel_properties:
        assert isinstance(getattr(channel, p), u.Quantity)
    # The preflight contamination is missing for the FUV channels and falls back to 1
    assert_quantity_allclose(channel._preflight_contamination, np.ones(channel.wavelength.shape))
