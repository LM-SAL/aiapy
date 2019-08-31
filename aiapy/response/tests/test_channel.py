"""
Tests for Channel class that holds response function and channel data
"""
import os
import collections

import pytest
import astropy.units as u
import astropy.time
from sunpy.util.metadata import MetaDict
try:
    import hissw
except ImportError:
    HISSW_MISSING = True
else:
    HISSW_MISSING = False

from aiapy.response import Channel
from aiapy.response.channel import VERSION_NUMBER

# FIXME: need a better way of testing whether to run SSW tests
needs_ssw = pytest.mark.skipif(
    HISSW_MISSING,
    reason='''hissw required for comparing SSW and Python results.
              IDL and SSW are also required.''')


@pytest.fixture
def channel():
    if HISSW_MISSING:
        instrument_file = None
    else:
        # Make sure we are using the same set of SSW data
        instrument_file = os.path.join(
            hissw.Environment().ssw_home,
            'sdo',
            'aia',
            'response',
            f'aia_V{VERSION_NUMBER}_all_fullinst.genx'
        )
    return Channel(94*u.angstrom, instrument_file=instrument_file)


@pytest.fixture
def hissw_env():
    return hissw.Environment(ssw_packages=['sdo/aia'], ssw_paths=['aia'])


def test_has_instrument_data(channel):
    assert hasattr(channel, '_instrument_data')
    assert isinstance(channel._instrument_data, collections.OrderedDict)


def test_has_channel_data(channel):
    assert hasattr(channel, '_data')
    assert isinstance(channel._data, MetaDict)


def test_channel_data_has_keys(channel):
    required_keys = ['wave', 'primary', 'secondary', 'fp_filter', 'ent_filter',
                     'geoarea', 'ccd', 'platescale', 'elecperev', 'elecperdn']
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


@pytest.mark.remote_data
def test_time_correction(channel):
    obstime = astropy.time.Time('2015-01-01T00:00:00', scale='utc')
    # NOTE: this just tests an expected result from aiapy, not necessarily an
    # absolutely correct result. It was calculated for the above time
    time_correction = 0.7667012041798814 * u.dimensionless_unscaled
    assert u.allclose(channel.time_correction(obstime), time_correction,
                      rtol=1e-10, atol=0.)


def test_eve_correction(channel):
    obstime = astropy.time.Time('2015-01-01T00:00:00', scale='utc')
    # NOTE: this just tests an expected result from aiapy, not necessarily an
    # absolutely correct result. It was calculated for the above time
    eve_correction = 1.0140518082508945 * u.dimensionless_unscaled
    assert u.allclose(channel.eve_correction(obstime), eve_correction,
                      rtol=1e-10, atol=0.)


def test_gain(channel):
    # FIXME: what is the best way to test this? do we even need to test this?
    ...


@needs_ssw
def test_wavelength_response_uncorrected(channel, hissw_env):
    r = channel.wavelength_response()
    ssw = hissw_env.run('r = aia_get_response(/area,/dn,evenorm=0)',
                        save_vars=['r'], verbose=False)
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)


@needs_ssw
def test_wavelength_response_no_crosstalk(channel, hissw_env):
    r = channel.wavelength_response(include_crosstalk=False)
    ssw = hissw_env.run('r = aia_get_response(/area,/dn,/noblend,evenorm=0)',
                        save_vars=['r'], verbose=False)
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)


@pytest.mark.remote_data
@needs_ssw
def test_wavelength_response_time(channel, hissw_env):
    now = astropy.time.Time.now()
    r = channel.wavelength_response(obstime=now)
    ssw = hissw_env.run(
        '''
r = aia_get_response(/area,/dn,evenorm=0,timedepend_date='{{obstime}}',version={{version}})
        ''',
        save_vars=['r'],
        args={'obstime': now.tai.isot, 'version': VERSION_NUMBER},
        verbose=False
    )
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)


@pytest.mark.remote_data
@needs_ssw
def test_wavelength_response_eve(channel, hissw_env):
    now = astropy.time.Time.now()
    r = channel.wavelength_response(obstime=now, include_eve_correction=True)
    ssw = hissw_env.run(
        '''
r = aia_get_response(/area,/dn,evenorm=1,timedepend_date='{{obstime}}',version={{version}})
        ''',
        save_vars=['r'],
        args={'obstime': now.tai.isot, 'version': VERSION_NUMBER},
        verbose=False
    )
    r_ssw = ssw['r'][f'A{channel.name}'][0]['ea'][0] * u.cm**2 * u.count / u.ph
    assert u.allclose(r, r_ssw, rtol=1e-4, atol=0. * u.cm**2 * u.count / u.ph)
