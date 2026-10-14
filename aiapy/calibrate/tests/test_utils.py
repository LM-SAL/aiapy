import re
import contextlib
import urllib.error
import urllib.request
from types import SimpleNamespace

import pytest

import astropy.units as u
from astropy.table import QTable
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.time import TimeRange

from aiapy.calibrate import degradation, utils
from aiapy.calibrate.utils import (
    _download_pointing_table,
    _select_epoch_from_correction_table,
    get_correction_table,
    get_error_table,
    get_pointing_table,
)
from aiapy.tests.data import get_test_filepath
from aiapy.utils.exceptions import AIApyUserWarning

# These are not fixtures so that they can be easily used in the parametrize mark
obstime = Time("2015-01-01T00:00:00", scale="utc")
correction_table_local = get_correction_table(get_test_filepath("aia_V8_20171210_050627_response_table.txt"))
CORRECTION_TABLE_COLUMNS = [
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


@pytest.fixture(scope="module")
def latest_jsoc_correction_table():
    return get_correction_table("jsoc")


@pytest.fixture(scope="module")
def latest_ssw_correction_table():
    return get_correction_table("ssw")


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
    for col in CORRECTION_TABLE_COLUMNS:
        assert col in table.colnames
    assert isinstance(table["T_START"], Time)
    assert isinstance(table["T_STOP"], Time)


@pytest.mark.parametrize(("wavelength", "wave_str"), [(94 * u.angstrom, "94_THIN"), (1600 * u.angstrom, "1600")])
def test_correction_table_selection(wavelength, wave_str) -> None:
    table = _select_epoch_from_correction_table(wavelength, obstime, correction_table_local, 8)
    assert isinstance(table, QTable)
    assert table.colnames == correction_table_local.colnames
    # The first epoch, then the epoch containing obstime, for one channel and version.
    assert len(table) == 2
    assert all(table["WAVE_STR"] == wave_str)
    assert all(table["VER_NUM"] == 8)
    assert table["T_START"][0] == table["T_START"].min()
    assert table["T_START"][1] <= obstime < table["T_STOP"][1]


def test_invalid_correction_table_input() -> None:
    with pytest.raises(
        ValueError,
        match=re.escape("correction_table must be a file path (pathlib.Path), 'JSOC' or 'SSW'. Not -1"),
    ):
        get_correction_table(source=-1)


def test_invalid_wavelength_raises_exception() -> None:
    with pytest.raises(ValueError, match=r'channel "1800.0 Angstrom" not in list of valid channels'):
        _select_epoch_from_correction_table(1800 * u.angstrom, obstime, correction_table_local)


def test_obstime_out_of_range() -> None:
    obstime_out_of_range = Time("2000-01-01T12:00:00", scale="utc")
    with pytest.raises(ValueError, match=f"No valid calibration epoch for {obstime_out_of_range}"):
        _select_epoch_from_correction_table(94 * u.angstrom, obstime_out_of_range, correction_table_local, 8)


@pytest.mark.remote_data
def test_pointing_table(pointing_table) -> None:
    expected_columns = ["T_START", "T_STOP"]
    for c in ["094", "131", "171", "193", "211", "304", "335", "1600", "1700", "4500"]:
        expected_columns += [
            f"A_{c}_X0",
            f"A_{c}_Y0",
            f"A_{c}_INSTROT",
            f"A_{c}_IMSCALE",
        ]
    t = Time("2011-01-01T00:00:00", scale="utc")
    table_jsoc = get_pointing_table("jsoc", time_range=(t - 3 * u.h, t + 3 * u.h))
    for table in [pointing_table, table_jsoc]:
        assert isinstance(table, QTable)
        for col in expected_columns:
            assert col in table.colnames
        assert isinstance(table["T_START"], Time)
        assert isinstance(table["T_STOP"], Time)
        # Ensure that none of the pointing parameters are masked columns
        for c in expected_columns[2:]:
            assert not hasattr(table[c], "mask")


@pytest.mark.parametrize(
    ("source", "time_range", "error", "match"),
    [
        ("jsoc", None, ValueError, "time_range must be provided when source is 'jsoc'"),
        ("jsoc", "2011-01-01", TypeError, "time_range must be a Time"),
        ("not a source", None, ValueError, "Invalid source: not a source"),
    ],
)
def test_pointing_table_invalid_input(source, time_range, error, match) -> None:
    with pytest.raises(error, match=match):
        get_pointing_table(source, time_range=time_range)


@pytest.mark.remote_data
def test_pointing_table_unavailable() -> None:
    # Check that missing pointing data raises a nice error
    t = Time("1990-01-01")
    with pytest.raises(RuntimeError, match="No data found for this query"):
        get_pointing_table("jsoc", time_range=Time([t - 3 * u.h, t + 3 * u.h]))


@pytest.mark.remote_data
def test_get_pointing_table_timerange() -> None:
    t = Time("2021-01-01T00:00:00", scale="utc")
    table = get_pointing_table("jsoc", time_range=TimeRange(t - 3 * u.h, t + 3 * u.h))
    assert isinstance(table, QTable)
    assert table["T_START"].min() == Time("2020-12-31T21:00:00", scale="utc")
    assert table["T_STOP"].max() == Time("2021-01-01T09:00:00", scale="utc")


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


@pytest.mark.remote_data
@pytest.mark.parametrize("channel", [94, 131, 171, 193, 211, 335] * u.angstrom)
def test_use_latest_version_in_correction_tables(
    channel, latest_jsoc_correction_table, latest_ssw_correction_table
) -> None:
    # Check that https://github.com/LM-SAL/aiapy/issues/375 is fixed
    # Issue is that we do not use the latest version of the correction table by default
    # if multiple versions are available for a given time
    time = Time("2010-11-03T12:15:00")
    d_jsoc = degradation(channel, time, correction_table=latest_jsoc_correction_table)
    d_ssw = degradation(channel, time, correction_table=latest_ssw_correction_table)
    assert_quantity_allclose(d_jsoc, d_ssw)


def _fake_urlopen(requests, *, body=b"", error=None):
    def urlopen(request, timeout):  # NOQA: ARG001
        requests.append(request)
        if error is not None:
            raise error
        return contextlib.nullcontext(SimpleNamespace(read=lambda: body, headers={"ETag": '"v2"'}))

    return urlopen


def test_download_pointing_table(monkeypatch, tmp_path) -> None:
    requests = []
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen(requests, body=b"new"))
    path = _download_pointing_table(tmp_path)
    assert path.read_bytes() == b"new"
    assert (tmp_path / "aia_pointing_table.etag").read_text() == '"v2"'
    assert requests[0].get_header("If-none-match") is None


def test_download_pointing_table_not_modified(monkeypatch, tmp_path) -> None:
    (tmp_path / "aia_pointing_table.csv").write_bytes(b"old")
    (tmp_path / "aia_pointing_table.etag").write_text('"v1"')
    requests = []
    not_modified = urllib.error.HTTPError(utils._POINTING_TABLE_URL, 304, "Not Modified", {}, None)
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen(requests, error=not_modified))
    assert _download_pointing_table(tmp_path).read_bytes() == b"old"
    assert requests[0].get_header("If-none-match") == '"v1"'


def test_download_pointing_table_offline(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(urllib.request, "urlopen", _fake_urlopen([], error=urllib.error.URLError("offline")))
    with pytest.raises(urllib.error.URLError, match="offline"):
        _download_pointing_table(tmp_path)
    (tmp_path / "aia_pointing_table.csv").write_bytes(b"old")
    with pytest.warns(AIApyUserWarning, match="Could not check for a newer pointing table"):
        assert _download_pointing_table(tmp_path).read_bytes() == b"old"


def test_lmsal_pointing_table_cached_per_session(monkeypatch, tmp_path) -> None:
    table_file = tmp_path / "aia_pointing_table.csv"
    table_file.write_text("T_START,T_STOP,A_094_X0\n2010-05-13T00:00:00Z,2010-05-13T03:00:00Z,2048.3\n")
    downloads = []

    def download():
        downloads.append(table_file)
        return table_file

    monkeypatch.setattr(utils, "_download_pointing_table", download)
    utils._lmsal_pointing_table.cache_clear()
    try:
        table = get_pointing_table("lmsal")
        table["A_094_X0"][0] = 0 * u.pix
        assert get_pointing_table("lmsal")["A_094_X0"][0] == 2048.3 * u.pix
        assert len(downloads) == 1
    finally:
        utils._lmsal_pointing_table.cache_clear()
