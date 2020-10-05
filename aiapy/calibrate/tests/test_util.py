"""
Tests for utility functions in response subpackage
"""
import pytest
import astropy.time
import astropy.table
import astropy.units as u

from aiapy.calibrate.util import get_correction_table, _select_epoch_from_table, get_pointing_table
from aiapy.tests.data import get_test_filepath


# These are not fixtures so that they can be easily used in the parametrize mark
obstime = astropy.time.Time('2015-01-01T00:00:00', scale='utc')
table_local = get_correction_table(
    correction_table=get_test_filepath('aia_V8_20171210_050627_response_table.txt'))


@pytest.mark.parametrize('correction_table', [
    pytest.param(None, marks=pytest.mark.remote_data),
    table_local,
    get_test_filepath('aia_V8_20171210_050627_response_table.txt'),
])
def test_correction_table(correction_table):
    table = get_correction_table(correction_table=correction_table)
    assert isinstance(table, astropy.table.QTable)
    expected_columns = ['VER_NUM',
                        'WAVE_STR',
                        'T_START',
                        'T_STOP',
                        'EFFA_P1',
                        'EFFA_P2',
                        'EFFA_P3',
                        'EFF_AREA',
                        'EFF_WVLN']
    assert all([cn in table.colnames for cn in expected_columns])
    assert isinstance(table['T_START'], astropy.time.Time)
    assert isinstance(table['T_STOP'], astropy.time.Time)


@pytest.mark.parametrize('wavelength', [94*u.angstrom, 1600*u.angstrom])
def test_correction_table_selection(wavelength):
    table = _select_epoch_from_table(wavelength,
                                     obstime,
                                     table_local,
                                     version=8)
    assert isinstance(table, astropy.table.QTable)
    expected_columns = ['VER_NUM',
                        'WAVE_STR',
                        'T_START',
                        'T_STOP',
                        'EFFA_P1',
                        'EFFA_P2',
                        'EFFA_P3',
                        'EFF_AREA',
                        'EFF_WVLN']
    assert all([cn in table.colnames for cn in expected_columns])
    assert isinstance(table['T_START'], astropy.time.Time)
    assert isinstance(table['T_STOP'], astropy.time.Time)


def test_invalid_correction_table_input():
    with pytest.raises(ValueError,
                       match='correction_table must be a file path, an existing table, or None.'):
        get_correction_table(correction_table=-1)


def test_invalid_wavelength_raises_exception():
    with pytest.raises(IndexError,
                       match='Correction table does not contain calibration for wavelength 1800'):
        _select_epoch_from_table(1800*u.angstrom, obstime, table_local)


def test_wrong_version_number_raises_exception():
    with pytest.raises(IndexError,
                       match='Correction table does not contain calibration for version -1'):
        _select_epoch_from_table(94*u.angstrom, obstime, table_local, version=-1)


def test_obstime_out_of_range():
    obstime_out_of_range = astropy.time.Time('2000-01-01T12:00:00', scale='utc')
    with pytest.raises(IndexError,
                       match=f'No valid calibration epoch for {obstime_out_of_range}'):
        _select_epoch_from_table(94*u.angstrom, obstime_out_of_range, table_local, version=8)


@pytest.mark.remote_data
def test_pointing_table():
    t = astropy.time.Time('2011-01-01T00:00:00', scale='utc')
    table = get_pointing_table(t-3*u.h, t+3*u.h)
    assert isinstance(table, astropy.table.QTable)
    expected_columns = ['T_START']
    for c in ['094', '171', '193', '211', '304', '335', '1600', '1700', '4500']:
        expected_columns += [f'A_{c}_X0', f'A_{c}_Y0', f'A_{c}_IMSCALE', f'A_{c}_IMSCALE']
    assert all([cn in table.colnames for cn in expected_columns])
    assert isinstance(table['T_START'], astropy.time.Time)
