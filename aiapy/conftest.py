import contextlib

import astropy.units as u
import pytest
import sunpy.data.test
import sunpy.map

ALL_CHANNELS = (94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500) * u.angstrom
CHANNELS = (94, 131, 171, 193, 211, 304, 335) * u.angstrom

# Force MPL to use non-gui backends for testing.
with contextlib.suppress(ImportError):
    import matplotlib

    matplotlib.use("Agg")


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
