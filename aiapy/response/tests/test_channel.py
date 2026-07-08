import collections
from pathlib import Path

import numpy as np
import pytest

import astropy.time
import astropy.units as u

from sunpy.util.metadata import MetaDict

from aiapy.calibrate.utils import get_correction_table
from aiapy.response import Channel
from aiapy.response.channel import VERSION_NUMBER
from aiapy.tests.data import get_test_filepath


# Mark all tests which use this fixture as online
@pytest.fixture(params=[pytest.param(None, marks=pytest.mark.remote_data)])
def channel(request, ssw_home):  # NOQA: ARG001
    if ssw_home is not None:
        instrument_file = Path(ssw_home) / "sdo" / "aia" / "response" / f"aia_V{VERSION_NUMBER}_all_fullinst.genx"
    else:
        instrument_file = None
    return Channel(94 * u.angstrom, instrument_file=instrument_file)


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


def test_has_instrument_data(channel) -> None:
    assert hasattr(channel, "_instrument_data")
    assert isinstance(channel._instrument_data, collections.OrderedDict)


def test_has_channel_data(channel) -> None:
    assert hasattr(channel, "_data")
    assert isinstance(channel._data, MetaDict)


def test_channel_data_has_keys(channel, required_keys) -> None:
    assert all(k in channel._data for k in required_keys)


def test_has_wavelength(channel) -> None:
    assert hasattr(channel, "wavelength")


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


def test_nominal_effective_area(channel):
    # Test that calculated effective area is the same as that stored in
    # the data file to within numerical precision
    channel.include_crosstalk = False
    assert u.allclose(
        channel.effective_area(),
        u.Quantity(channel._data["effarea"], "cm2"),
        atol=None,
        rtol=1e-6,
    )


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
    correction_table = get_correction_table(source=source)
    channel.correction_table = correction_table
    channel.calibration_version = np.max(correction_table["VER_NUM"])
    obstime = astropy.time.Time("2015-01-01T00:00:00", scale="utc")
    eve_correction = channel._get_eve_correction(obstime)
    assert u.allclose(eve_correction, eve_correction_truth, rtol=1e-6, atol=None)


@pytest.mark.parametrize(
    ("obstime", "correction_table", "calibration_version", "include_eve_correction", "include_crosstalk"),
    [
        (None, None, None, False, True),
        (None, None, None, False, False),
        ("2020-01-01", get_test_filepath("aia_V8_20171210_050627_response_table.txt"), 8, False, True),
        (astropy.time.Time.now(), get_test_filepath("aia_V8_20171210_050627_response_table.txt"), 8, True, True),
    ],
)
def test_effective_area(
    channel, obstime, correction_table, calibration_version, include_eve_correction, include_crosstalk
):
    # NOTE: this does not test correctness, but just that the method can
    # be run with the various combinations of inputs. The tests below test the
    # correctness of the output as evaluated by their similarity to those
    # results from SSW.
    channel.include_crosstalk = include_crosstalk
    channel.include_eve_correction = include_eve_correction
    channel.calibration_version = calibration_version
    channel.correction_table = correction_table
    bound = channel if obstime is None else channel.at(obstime)
    assert bound.effective_area().shape == channel.wavelength.shape


def test_wavelength_response_uncorrected(channel, idl_environment) -> None:
    channel.correction_table = get_correction_table(
        source=get_test_filepath("aia_V8_20171210_050627_response_table.txt")
    )
    r = channel.wavelength_response()
    ssw = idl_environment.run("r = aia_get_response(/area,/dn,evenorm=0)", save_vars=["r"], verbose=False)
    r_ssw = ssw["r"][f"A{channel.name}"][0]["ea"][0] * u.cm**2 * u.DN / u.ph
    r_ssw *= channel.pixel_solid_angle
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=None)


def test_wavelength_response_no_crosstalk(channel, idl_environment):
    channel.include_crosstalk = False
    r = channel.wavelength_response()
    ssw = idl_environment.run(
        "r = aia_get_response(/area,/dn,/noblend,evenorm=0)",
        save_vars=["r"],
        verbose=False,
    )
    r_ssw = ssw["r"][f"A{channel.name}"][0]["ea"][0] * u.cm**2 * u.DN / u.ph
    r_ssw *= channel.pixel_solid_angle
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=None)


@pytest.mark.parametrize("include_eve_correction", [False, True])
def test_wavelength_response_time(channel, idl_environment, include_eve_correction) -> None:
    now = astropy.time.Time.now()
    correction_table = get_test_filepath("aia_V8_20171210_050627_response_table.txt")
    calibration_version = 8
    channel.calibration_version = calibration_version
    channel.correction_table = correction_table
    channel.include_eve_correction = include_eve_correction
    channel.include_crosstalk = True
    r = channel.at(now).wavelength_response()
    ssw = idl_environment.run(
        """
        r = aia_get_response(/area,/dn,evenorm={{evenorm}}, $
        timedepend_date='{{obstime}}',version={{version}},respversion='{{respversion}}')
        """,
        save_vars=["r"],
        args={
            "obstime": now.tai.isot,
            "evenorm": int(include_eve_correction),
            "version": 8,
            "respversion": "8",
        },
        verbose=False,
    )
    r_ssw = ssw["r"][f"A{channel.name}"][0]["ea"][0] * u.cm**2 * u.DN / u.ph
    r_ssw *= channel.pixel_solid_angle
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=None)


@pytest.mark.remote_data
@pytest.mark.parametrize("channel_wavelength", [1600 * u.angstrom, 1700 * u.angstrom, 4500 * u.angstrom])
def test_fuv_channel(channel_wavelength, channel_properties, required_keys) -> None:
    # There are a few corner cases for the 1600, 1700, and 4500 channels
    channel = Channel(channel_wavelength, include_eve_correction=False)
    assert all(k in channel._data for k in required_keys)
    for p in channel_properties:
        assert isinstance(getattr(channel, p), u.Quantity)
    assert u.allclose(channel.degradation(), 1, rtol=0.0, atol=None)
