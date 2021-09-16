import pytest

import astropy.units as u
import sunpy.data.test
import sunpy.map

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

# Do not require hissw for tests
try:
    import hissw
except ImportError:
    pass


@pytest.fixture
def aia_171_map():
    m = sunpy.map.Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits'))
    # For testing purposes, need the map to be 4K-by-4K
    return m.resample((4096, 4096)*u.pixel)


@pytest.fixture(scope='session')
def idl_environment():
    if idl_available():
        return hissw.Environment(ssw_packages=['sdo/aia'], ssw_paths=['aia'])
    else:
        pytest.skip(
            "A working IDL installation is not available. You will not be able to run portions of the test suite."
        )


@pytest.fixture(scope='session')
def ssw_home():
    if idl_available():
        return hissw.Environment().ssw_home
    return None


def idl_available():
    try:
        import hissw
        _ = hissw.Environment().run('')
        return True
    except Exception:
        return False
