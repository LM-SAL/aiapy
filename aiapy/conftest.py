import matplotlib as mpl

# Force MPL to use non-gui backends for testing.
mpl.use("Agg")

import numpy as np
import pytest

from astropy.io import fits

import sunpy.data.test
import sunpy.map
from sunpy import log

from aiapy.tests.helpers import CHANNELS
from aiapy.utils import detector_dimensions


@pytest.fixture
def aia_171_test_map():
    """The real AIA level-1 observation distributed by SunPy, resampled to 128 pixels."""
    return sunpy.map.Map(sunpy.data.test.get_test_filepath("aia_171_level1.fits"))


@pytest.fixture
def aia_171_map(aia_171_test_map):
    """Upsample the test observation for algorithms requiring the full detector shape."""
    return aia_171_test_map.resample(detector_dimensions())


@pytest.fixture(params=["hmi_sharp_magnetogram.header", "hmi_bharp_vlos_mag.header"])
def hmi_map(request):
    """Synthetic pixels with the header and geometry of real HMI cutout/full-frame files."""
    header = fits.Header.fromtextfile(sunpy.data.test.get_test_filepath(request.param))
    # Floating control pixels do not use the source file's integer blank sentinel.
    header.remove("BLANK", ignore_missing=True)
    data = np.ones((header["naxis2"], header["naxis1"]), dtype=np.float32)
    hdu = fits.PrimaryHDU(data, header=header)
    return sunpy.map.Map(hdu.data, hdu.header)


@pytest.fixture
def aia_193_level1_map():
    # Need an actual 4K-by-4K map to do the spike replacement
    return sunpy.map.Map(
        "https://github.com/sunpy/data/raw/refs/heads/main/aiapy/aia_lev1_193a_2013_03_15t12_01_06_84z_image_lev1.fits",
    )


@pytest.fixture
def psf(aia_171_test_map):
    import aiapy.psf  # NOQA: PLC0415

    return aiapy.psf.calculate_psf(aia_171_test_map.wavelength, use_preflightcore=True, diffraction_orders=[-1, 0, 1])


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
