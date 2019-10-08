"""
Tests for utility functions in response subpackage
"""
import pytest
import astropy.time
import astropy.table

from aiapy.calibrate.util import get_correction_table
from aiapy.tests.data import get_test_filepath


@pytest.mark.remote_data
@pytest.fixture
def correction_table_jsoc():
    return get_correction_table()


@pytest.fixture
def correction_table_file():
    return get_correction_table(
        correction_table=get_test_filepath(
            'aia_V8_20171210_050627_response_table.txt'
        ))


@pytest.fixture
def expected_columns():
    return [
        'VER_NUM',
        'WAVE_STR',
        'T_START',
        'T_STOP',
        'EFFA_P1',
        'EFFA_P2',
        'EFFA_P3',
        'EFF_AREA',
        'EFF_WVLN'
    ]


@pytest.mark.remote_data
def test_correction_table_jsoc_is_table(correction_table_jsoc):
    assert isinstance(correction_table_jsoc, astropy.table.Table)


def test_correction_table_file_is_table(correction_table_file):
    assert isinstance(correction_table_file, astropy.table.Table)


@pytest.mark.remote_data
def test_correction_table_jsoc_has_columns(correction_table_jsoc,
                                           expected_columns):
    assert all([cn in correction_table_jsoc.colnames
                for cn in expected_columns])


def test_correction_table_file_has_columns(correction_table_file,
                                           expected_columns):
    assert all([cn in correction_table_file.colnames
                for cn in expected_columns])


@pytest.mark.remote_data
def test_correction_table_jsoc_has_time_columns(correction_table_jsoc):
    assert isinstance(correction_table_jsoc['T_START'], astropy.time.Time)
    assert isinstance(correction_table_jsoc['T_STOP'], astropy.time.Time)


def test_correction_table_file_has_time_columns(correction_table_file):
    assert isinstance(correction_table_file['T_START'], astropy.time.Time)
    assert isinstance(correction_table_file['T_STOP'], astropy.time.Time)
