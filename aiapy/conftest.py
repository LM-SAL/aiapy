import contextlib

import pytest
from numpy.random import default_rng

import astropy.units as u
from astropy.io import fits

import sunpy.data.test
import sunpy.map
from sunpy import log

RANDOM_GENERATOR = default_rng()
CHANNELS = [94, 131, 171, 193, 211, 304, 335] * u.angstrom
ALL_CHANNELS = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500] * u.angstrom

with contextlib.suppress(ImportError):
    import matplotlib as mpl

    # Force MPL to use non-gui backends for testing.
    mpl.use("Agg")


@pytest.fixture
def aia_171_map():
    m = sunpy.map.Map(sunpy.data.test.get_test_filepath("aia_171_level1.fits"))
    # For testing purposes, need the map to be 4K-by-4K
    return m.resample((4096, 4096) * u.pixel)


@pytest.fixture
def aia_193_level1_map():
    # Need an actual 4K-by-4K map to do the spike replacement
    return sunpy.map.Map(
        "https://github.com/sunpy/data/blob/main/aiapy/aia_lev1_193a_2013_03_15t12_01_06_84z_image_lev1.fits?raw=true",
    )


@pytest.fixture
def psf_94():
    import aiapy.psf  # NOQA: PLC0415

    return aiapy.psf.calculate_psf(CHANNELS[0], use_preflightcore=True)


@pytest.fixture
def remote_131_map():
    return sunpy.map.Map(
        "https://github.com/sunpy/data/raw/refs/heads/main/aiapy/aia.lev1_euv_12s.2024-05-08T014308Z.131.image_lev1.fits"
    )


@pytest.fixture
def remote_131_idl_deconvolved_output():
    # IDL was run with: aia_deconvolve_richardsonlucy(image,psf,niter=30)
    return fits.getdata("https://github.com/sunpy/data/raw/refs/heads/main/aiapy/deconvolved_image_131.fits.fz")


def idl_available() -> bool | None:
    try:
        import hissw  # NOQA: PLC0415

        hissw.Environment().run("")
        return True  # NOQA: TRY300
    except Exception as e:  # NOQA: BLE001
        log.warning(e)
        return False


@pytest.fixture(scope="session")
def idl_environment():
    if idl_available():
        import hissw  # NOQA: PLC0415

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
