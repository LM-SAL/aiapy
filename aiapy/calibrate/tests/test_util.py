import re

import pytest

import astropy.units as u
from astropy.table import QTable
from astropy.time import Time

from aiapy.calibrate.util import (
    _select_epoch_from_correction_table,
    get_correction_table,
    get_error_table,
    get_pointing_table,
)
from aiapy.tests.data import get_test_filepath

# These are not fixtures so that they can be easily used in the parametrize mark
obstime = Time("2015-01-01T00:00:00", scale="utc")
correction_table_local = get_correction_table(get_test_filepath("aia_V8_20171210_050627_response_table.txt"))


@pytest.mark.parametrize(
    "source",
    [
        pytest.param("JSoC", marks=pytest.mark.remote_data),  # To check the lower case comparison is working
        pytest.param("SsW", marks=pytest.mark.remote_data),  # To check the lower case comparison is workings
        get_test_filepath("aia_V8_20171210_050627_response_table.txt"),
    ],
)
def test_correction_table(source) -> None:
    table = get_correction_table(source=source)
    assert isinstance(table, QTable)
    expected_columns = [
        "VER_NUM",
        "WAVE_STR",
        "T_START",
        "T_STOP",
        "EFFA_P1",
        "EFFA_P2",
        "EFFA_P3",
        "EFF_AREA",
        "EFF_WVLN",
    ]
    assert all(cn in table.colnames for cn in expected_columns)
    assert isinstance(table["T_START"], Time)
    assert isinstance(table["T_STOP"], Time)


@pytest.mark.parametrize("wavelength", [94 * u.angstrom, 1600 * u.angstrom])
def test_correction_table_selection(wavelength) -> None:
    table = _select_epoch_from_correction_table(wavelength, obstime, correction_table_local)
    assert isinstance(table, QTable)
    expected_columns = [
        "VER_NUM",
        "WAVE_STR",
        "T_START",
        "T_STOP",
        "EFFA_P1",
        "EFFA_P2",
        "EFFA_P3",
        "EFF_AREA",
        "EFF_WVLN",
    ]
    assert all(cn in table.colnames for cn in expected_columns)
    assert isinstance(table["T_START"], Time)
    assert isinstance(table["T_STOP"], Time)


def test_invalid_correction_table_input() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("correction_table must be a file path (pathlib.Path), 'JSOC' or 'SSW'. Not -1"),
    ):
        get_correction_table(source=-1)


def test_invalid_wavelength_raises_exception() -> None:
    with pytest.raises(ValueError, match='channel "1800.0 Angstrom" not in list of valid channels'):
        _select_epoch_from_correction_table(1800 * u.angstrom, obstime, correction_table_local)


def test_obstime_out_of_range() -> None:
    obstime_out_of_range = Time("2000-01-01T12:00:00", scale="utc")
    with pytest.raises(ValueError, match=f"No valid calibration epoch for {obstime_out_of_range}"):
        _select_epoch_from_correction_table(94 * u.angstrom, obstime_out_of_range, correction_table_local)


@pytest.mark.remote_data
def test_pointing_table() -> None:
    expected_columns = ["T_START", "T_STOP"]
    for c in ["094", "171", "193", "211", "304", "335", "1600", "1700", "4500"]:
        expected_columns += [
            f"A_{c}_X0",
            f"A_{c}_Y0",
            f"A_{c}_IMSCALE",
            f"A_{c}_IMSCALE",
        ]
    t = Time("2011-01-01T00:00:00", scale="utc")
    table_lmsal = get_pointing_table("lmsal")
    table_jsoc = get_pointing_table("jsoc", time_range=(t - 3 * u.h, t + 3 * u.h))
    for table in [table_lmsal, table_jsoc]:
        assert isinstance(table, QTable)
        assert all(cn in table.colnames for cn in expected_columns)
        assert isinstance(table["T_START"], Time)
        assert isinstance(table["T_STOP"], Time)
        # Ensure that none of the pointing parameters are masked columns
        for c in expected_columns[2:]:
            assert not hasattr(table[c], "mask")


@pytest.mark.remote_data
def test_pointing_table_unavailable() -> None:
    # Check that missing pointing data raises a nice error
    t = Time("1990-01-01")
    with pytest.raises(RuntimeError, match="No data found for this query"):
        get_pointing_table("jsoc", time_range=Time([t - 3 * u.h, t + 3 * u.h]))


@pytest.mark.parametrize(
    "error_table",
    [
        pytest.param("SSW", marks=pytest.mark.remote_data),
        get_test_filepath("aia_V3_error_table.txt"),
        str(get_test_filepath("aia_V3_error_table.txt")),
    ],
)
def test_error_table(error_table) -> None:
    table = get_error_table(error_table)
    assert isinstance(table, QTable)
    assert len(table) == 10


def test_invalid_error_table_input() -> None:
    with pytest.raises(TypeError, match="source must be a filepath, or 'SSW', not -1"):
        get_error_table(-1)
