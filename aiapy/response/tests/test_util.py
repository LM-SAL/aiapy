"""
Tests for utility functions in response subpackage
"""
import pytest
import astropy.time
import astropy.table

from aiapy.response.util import get_correction_table


@pytest.mark.remote_data
@pytest.fixture
def correction_table():
    return get_correction_table()


@pytest.mark.remote_data
def test_correction_table_is_table(correction_table):
    assert isinstance(correction_table, astropy.table.Table)


@pytest.mark.remote_data
def test_correction_table_has_columns(correction_table):
    colnames = ['VER_NUM', 'WAVE_STR', 'T_START', 'T_STOP', 'EFFA_P1', 'EFFA_P2', 'EFFA_P3',
                'EFF_AREA', 'EFF_WVLN']
    assert all([cn in correction_table.colnames for cn in colnames])


@pytest.mark.remote_data
def test_correction_table_has_time_columns(correction_table):
    assert isinstance(correction_table['T_START'], astropy.time.Time)
    assert isinstance(correction_table['T_STOP'], astropy.time.Time)
