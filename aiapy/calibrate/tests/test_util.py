"""
Tests for utility functions in response subpackage
"""
import pytest
import astropy.time
import astropy.table
import astropy.units as u

from aiapy.calibrate.util import get_correction_table, _select_epoch_from_table
from aiapy.tests.data import get_test_filepath


# These are not fixtures so that they can be easily used in the parametrize mark
obstime = astropy.time.Time('2015-01-01T00:00:00')
table_local = get_correction_table(
    correction_table=get_test_filepath('aia_V8_20171210_050627_response_table.txt'))


@pytest.mark.parametrize('correction_table', [
    pytest.param(None, marks=pytest.mark.remote_data),
    table_local,
    _select_epoch_from_table(94*u.angstrom,
                             obstime,
                             correction_table=table_local,
                             calibration_version=8),
    _select_epoch_from_table(1600*u.angstrom,
                             obstime,
                             correction_table=table_local,
                             calibration_version=8),
])
def test_correction_table(correction_table):
    # NOTE: This is kind of a hack to get around the fact that sticking the table
    # generated from the remote data will get triggered regardless of the remote
    # mark.
    if correction_table is None:
        correction_table = get_correction_table()
    assert isinstance(correction_table, astropy.table.Table)
    expected_columns = ['VER_NUM',
                        'WAVE_STR',
                        'T_START',
                        'T_STOP',
                        'EFFA_P1',
                        'EFFA_P2',
                        'EFFA_P3',
                        'EFF_AREA',
                        'EFF_WVLN']
    assert all([cn in correction_table.colnames for cn in expected_columns])
    assert isinstance(correction_table['T_START'], astropy.time.Time)
    assert isinstance(correction_table['T_STOP'], astropy.time.Time)


def test_invalid_wavelength_raises_exception():
    with pytest.raises(IndexError,
                       match='Correction table does not contain calibration for wavelength 1800'):
        _select_epoch_from_table(1800*u.angstrom,
                                 obstime,
                                 correction_table=table_local)


def test_wrong_version_number_raises_exception():
    with pytest.raises(IndexError,
                       match='Correction table does not contain calibration for version -1'):
        _select_epoch_from_table(94*u.angstrom,
                                 obstime,
                                 correction_table=table_local,
                                 calibration_version=-1)
