"""
Tests for Channel class that holds response function and channel data
"""
import os
import collections

import pytest
import astropy.units as u
import astropy.time
from sunpy.util.metadata import MetaDict

from aiapy.response import Channel
from aiapy.response.channel import VERSION_NUMBER
from aiapy.calibrate.util import get_correction_table
from aiapy.tests.data import get_test_filepath


@pytest.fixture
def channel(ssw_home):
    if ssw_home is not None:
        instrument_file = os.path.join(
            ssw_home,
            'sdo',
            'aia',
            'response',
            f'aia_V{VERSION_NUMBER}_all_fullinst.genx'
        )
    else:
        instrument_file = None

    return Channel(94*u.angstrom, instrument_file=instrument_file)


def test_has_instrument_data(channel):
    assert hasattr(channel, '_instrument_data')
    assert isinstance(channel._instrument_data, collections.OrderedDict)


def test_has_channel_data(channel):
    assert hasattr(channel, '_data')
    assert isinstance(channel._data, MetaDict)


def test_channel_data_has_keys(channel):
    required_keys = [
        'wave',
        'primary',
        'secondary',
        'fp_filter',
        'ent_filter',
        'geoarea',
        'ccd',
        'platescale',
        'elecperev',
        'elecperdn',
    ]
    assert all([k in channel._data for k in required_keys])


def test_has_wavelength(channel):
    assert hasattr(channel, 'wavelength')


def test_channel_wavelength(channel):
    assert channel.channel == 94*u.angstrom
    assert channel.name == '94'


def test_channel_properties(channel):
    """Test that expected properties are present and are quantities"""
    properties = [
        'wavelength',
        'primary_reflectance',
        'secondary_reflectance',
        'focal_plane_filter_efficiency',
        'entrance_filter_efficiency',
        'geometrical_collecting_area',
        'quantum_efficiency',
        'contamination',
        'plate_scale',
        'effective_area',
        'gain'
    ]
    for p in properties:
        assert isinstance(getattr(channel, p), u.Quantity)


def test_effective_area(channel):
    effective_area = (channel.primary_reflectance
                      * channel.secondary_reflectance
                      * channel.focal_plane_filter_efficiency
                      * channel.entrance_filter_efficiency
                      * channel.geometrical_collecting_area
                      * channel.quantum_efficiency
                      * channel.contamination)
    assert (effective_area == channel.effective_area).all()


@pytest.mark.parametrize('correction_table,eve_correction_truth', [
    pytest.param(None, 1.0140518082508945 * u.dimensionless_unscaled,
                 marks=pytest.mark.remote_data),
    (get_test_filepath('aia_V8_20171210_050627_response_table.txt'),
     1.0140386988603103 * u.dimensionless_unscaled),
    (get_correction_table(correction_table=get_test_filepath(
        'aia_V8_20171210_050627_response_table.txt')),
     1.0140386988603103 * u.dimensionless_unscaled),
])
def test_eve_correction(channel, correction_table, eve_correction_truth):
    # NOTE: this just tests an expected result from aiapy, not necessarily an
    # absolutely correct result. It was calculated for the above time and
    # the correction parameters in JSOC at the time this code was committed/
    # the specific correction table file.
    # NOTE: If the first test starts failing, it may be because the correction
    # table parameters have been updated in JSOC.
    # NOTE: This second and third values are different from first because the
    # first is the expected result from JSOC and the second and third are
    # results from the correction table file in SSW. The result returned by
    # JSOC are not necessarily the same as those in the correction table files
    # in SSW though they should be close.
    obstime = astropy.time.Time('2015-01-01T00:00:00', scale='utc')
    eve_correction = channel.eve_correction(
        obstime, correction_table=correction_table)
    assert u.allclose(eve_correction, eve_correction_truth,
                      rtol=1e-10, atol=0.)


def test_gain(channel):
    # FIXME: what is the best way to test this? do we even need to test this?
    ...


def test_wavelength_response_uncorrected(channel, idl_environment):
    r = channel.wavelength_response()
    ssw = idl_environment.run('r = aia_get_response(/area,/dn,evenorm=0)',
                              save_vars=['r'],
                              verbose=False)
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)


def test_wavelength_response_no_crosstalk(channel, idl_environment):
    r = channel.wavelength_response(include_crosstalk=False)
    ssw = idl_environment.run(
        'r = aia_get_response(/area,/dn,/noblend,evenorm=0)',
        save_vars=['r'],
        verbose=False
    )
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)


@pytest.mark.parametrize('include_eve_correction', [False, True])
def test_wavelength_response_time(channel, idl_environment,
                                  include_eve_correction):
    now = astropy.time.Time.now()
    correction_table = get_test_filepath(
        'aia_V8_20171210_050627_response_table.txt')
    r = channel.wavelength_response(
        obstime=now,
        include_eve_correction=include_eve_correction,
        correction_table=correction_table)
    ssw = idl_environment.run(
        '''
        r = aia_get_response(/area,/dn,evenorm={{evenorm}}, $
        timedepend_date='{{obstime}}',version={{version}},respversion='{{respversion}}')
        ''',
        save_vars=['r'],
        args={'obstime': now.tai.isot,
              'version': VERSION_NUMBER,
              'evenorm': int(include_eve_correction),
              'respversion': '20171210_050627'},
        verbose=False
    )
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)
