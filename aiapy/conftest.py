import pytest
import hissw
from astropy.version import version as astropy_version
if astropy_version < '3.0':
    # For older versions of astropy
    from astropy.tests.pytest_plugins import *
else:
    from astropy.tests.plugins.display import (PYTEST_HEADER_MODULES,
                                               TESTED_VERSIONS)
from astropy.tests.helper import enable_deprecations_as_exceptions


@pytest.fixture(scope='session')
def idl_environment():
    if idl_available():
        return hissw.Environment(ssw_packages=['sdo/aia'], ssw_paths=['aia'])
    else:
        pytest.skip('''A working IDL installation is not available. You will
                       not be able to run portions of the test suite.''')


@pytest.fixture(scope='session')
def ssw_home():
    if idl_available():
        return hissw.Environment().ssw_home
    else:
        return None


def idl_available():
    try:
        _ = hissw.Environment().run('')
    except (FileNotFoundError, ValueError):
        return False
    else:
        return True
