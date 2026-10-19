from pathlib import Path

import matplotlib as mpl

# Force MPL to use non-gui backends for testing.
mpl.use("Agg")

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

import sunpy.data.test
import sunpy.map
from sunpy import log
from sunpy.data import cache

from aiapy.calibrate.utils import get_correction_table, get_pointing_table
from aiapy.tests.helpers import CHANNELS, make_aia_map
from aiapy.utils import detector_dimensions


@pytest.fixture
def aia_171_test_map():
    """
    SunPy's 128x128 copy of a real AIA level-1 observation.
    """
    return sunpy.map.Map(sunpy.data.test.get_test_filepath("aia_171_level1.fits"))


@pytest.fixture
def aia_171_map(aia_171_test_map):
    """
    Upsample the test observation for algorithms requiring the full detector shape.
    """
    return aia_171_test_map.resample(detector_dimensions())


@pytest.fixture
def aia_171_submap(aia_171_map):
    top_right = SkyCoord(0 * u.arcsec, -200 * u.arcsec, frame=aia_171_map.coordinate_frame)
    bottom_left = SkyCoord(-900 * u.arcsec, -900 * u.arcsec, frame=aia_171_map.coordinate_frame)
    return aia_171_map.submap(bottom_left, top_right=top_right)


@pytest.fixture
def small_aia_map():
    coordinate = SkyCoord(
        600 * u.arcsec, -200 * u.arcsec, frame="helioprojective", observer="earth", obstime="2020-01-01"
    )
    return make_aia_map(
        np.arange(32 * 40, dtype=np.int16).reshape(32, 40),
        171,
        unit=u.DN,
        level=1,
        rotation=23.0,
        coordinate=coordinate,
    )


@pytest.fixture(params=["hmi_sharp_magnetogram.header", "hmi_bharp_vlos_mag.header"])
def hmi_map(request):
    """
    Synthetic pixels with the header and geometry of real HMI cutout/full-frame files.
    """
    header = fits.Header.fromtextfile(sunpy.data.test.get_test_filepath(request.param))
    # Floating control pixels do not use the source file's integer blank sentinel.
    header.remove("BLANK", ignore_missing=True)
    data = np.ones((header["naxis2"], header["naxis1"]), dtype=np.float32)
    hdu = fits.PrimaryHDU(data, header=header)
    return sunpy.map.Map(hdu.data, hdu.header)


@pytest.fixture(scope="session")
def aia_193_level1_file():
    return cache.download(
        "https://github.com/sunpy/data/raw/refs/heads/main/aiapy/aia_lev1_193a_2013_03_15t12_01_06_84z_image_lev1.fits",
    )


@pytest.fixture(scope="session")
def aia_193_level1_map(aia_193_level1_file):
    # Need an actual 4K-by-4K map to do the spike replacement.
    # Shared between tests, so do not modify it in a test.
    return sunpy.map.Map(aia_193_level1_file)


@pytest.fixture(scope="session")
def jsoc_correction_table():
    # Shared between tests so the JSOC is queried once
    return get_correction_table("jsoc")


@pytest.fixture(scope="session")
def pointing_table():
    # Shared between tests, so copy it before modifying it in a test.
    return get_pointing_table("lmsal")


@pytest.fixture
def psf(aia_171_test_map):
    import aiapy.psf  # NOQA: PLC0415

    return aiapy.psf.calculate_psf(aia_171_test_map.wavelength, use_preflightcore=True, diffraction_orders=[-1, 0, 1])


@pytest.fixture
def psf_128(psf):
    """
    The center of ``psf``, renormalized to match the 128x128 test map.
    """
    center = psf[2048 - 64 : 2048 + 64, 2048 - 64 : 2048 + 64]
    return center / center.sum()


@pytest.fixture
def aia_171_positive_map(aia_171_test_map):
    """
    The test map clipped to at least 1 DN.

    SSW's Richardson-Lucy keeps negative pixels and returns the absolute value, while
    aiapy clips them to zero, so the two are only comparable for positive images.
    """
    data = np.clip(aia_171_test_map.data, 1, None).astype(float)
    return aia_171_test_map._new_instance(data, aia_171_test_map.meta)


@pytest.fixture
def psf_94():
    import aiapy.psf  # NOQA: PLC0415

    return aiapy.psf.calculate_psf(CHANNELS[0], use_preflightcore=True)


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
            ssw_packages=["sdo/aia", "vobs/ontology"],
            ssw_paths=["aia", "ontology"],
            extra_paths=[str(Path(__file__).parent / "tests" / "idl")],
        )
    pytest.skip(
        "A working IDL installation is not available. You will not be able to run portions of the test suite.",
    )
