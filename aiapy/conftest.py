import contextlib

import astropy.units as u
import pytest
import sunpy.data.test
import sunpy.map
from sunpy import log

ALL_CHANNELS = (94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500) * u.angstrom
CHANNELS = (94, 131, 171, 193, 211, 304, 335) * u.angstrom

# Force MPL to use non-gui backends for testing.
with contextlib.suppress(ImportError):
    import matplotlib as mpl

    mpl.use("Agg")


@pytest.fixture()
def aia_171_map():
    m = sunpy.map.Map(sunpy.data.test.get_test_filepath("aia_171_level1.fits"))
    # For testing purposes, need the map to be 4K-by-4K
    return m.resample((4096, 4096) * u.pixel)


@pytest.fixture()
def all_channels():
    return CHANNELS


@pytest.fixture()
def channels():
    return CHANNELS


@pytest.fixture()
def psf_94(channels):
    import aiapy.psf

    return aiapy.psf.psf(channels[0], use_preflightcore=True)


def idl_available():
    try:
        import hissw

        hissw.Environment().run("")
        return True  # NOQA: TRY300
    except Exception as e:  # NOQA: BLE001
        log.warning(e)
        return False


@pytest.fixture(scope="session")
def idl_environment():
    if idl_available():
        import hissw

        return hissw.Environment(
            ssw_packages=["sdo/aia"],
            ssw_paths=["aia"],
        )
    pytest.skip(
        "A working IDL installation is not available. You will not be able to run portions of the test suite.",
    )


@pytest.fixture(scope="session")
def ssw_home(idl_environment):
    return idl_environment.ssw_home if idl_available() else None
