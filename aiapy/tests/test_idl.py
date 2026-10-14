"""
Contains all the the IDL specific tests for aiapy.
"""

from pathlib import Path

import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.tests.helper import assert_quantity_allclose

import aiapy.psf
from aiapy.calibrate import estimate_error, register
from aiapy.calibrate.utils import get_error_table
from aiapy.tests.helpers import CHANNELS


@pytest.mark.parametrize(
    ("channel", "counts", "include_eve", "include_preflight", "include_chianti"),
    [[c, 10 * u.DN / u.pixel] + 3 * [False] for c in CHANNELS]
    + [
        [171 * u.angstrom, 1000 * u.DN / u.pix, True, False, False],
        [171 * u.angstrom, 1000 * u.DN / u.pix, False, False, True],
    ],
)
def test_error_consistent(idl_environment, channel, counts, include_eve, include_preflight, include_chianti) -> None:
    idl = """
    common aia_bp_error_common,common_errtable
    common_errtable=aia_bp_read_error_table('{{ error_table }}')
    data = {{ data }}
    channel = {{ channel }}
    error=aia_bp_estimate_error(data,channel,n_sample=1{{ include_eve }}{{ include_preflight }}{{ include_chianti }})
    """
    error_table_path = Path(idl_environment.ssw_home) / "sdo" / "aia" / "response" / "aia_V3_error_table.txt"
    error_table = get_error_table(error_table_path)
    ssw = idl_environment.run(
        idl,
        save_vars=["error"],
        args={
            "channel": channel.to("angstrom").value,
            "data": counts.to("dn pixel-1").value,
            "error_table": error_table_path.as_posix(),
            "include_eve": ",/evenorm" if include_eve else "",
            # NOTE: use of this keyword is actually broken in SSW so these
            # tests only set it to False until it works, but consistency with
            # these results has been verified
            "include_preflight": ",/cal" if include_preflight else "",
            "include_chianti": ",/temperature" if include_chianti else "",
        },
        verbose=False,
    )
    error_ssw = ssw["error"] * counts.unit
    error = estimate_error(
        counts,
        channel,
        include_eve=include_eve,
        include_preflight=include_preflight,
        include_chianti=include_chianti,
        error_table=error_table,
        compare_idl=True,
    )
    assert_quantity_allclose(error, error_ssw, rtol=1e-4)


# For now, we only need the 94 Angstrom PSF
# In future, we will want to test all channels
@pytest.fixture(params=[CHANNELS[0]])
def psf_idl(idl_environment, request):
    """
    The point spread function as calculated by aia_calc_psf.pro.
    """
    r = idl_environment.run(
        "psf = aia_calc_psf({{channel}},/use_preflightcore)",
        args={"channel": f"{request.param.value:.0f}"},
        save_vars=["psf"],
        verbose=False,
    )
    return r["psf"]


def test_psf_consistent(psf_94, psf_idl) -> None:
    """
    Check whether PSF is consistent with IDL calculation.

    **WARNING** This test will take a very long time to run.
    """
    # NOTE: Only compare values above some threshold as the
    # rest of the PSF is essentially noise
    i_valid = np.where(psf_idl > 1e-10)
    # NOTE: The IDL and Python PSF functions have been found to
    # agree within 0.2% for all points along the PSF arms for
    # both the preflight and non-preflight cases.
    np.testing.assert_allclose(psf_94[i_valid], psf_idl[i_valid], rtol=2e-3)


def test_deconvolve_consistent(idl_environment, aia_171_positive_map, psf_128, tmp_path) -> None:
    fits.writeto(tmp_path / "image.fits", aia_171_positive_map.data)
    fits.writeto(tmp_path / "psf.fits", psf_128)
    ssw = idl_environment.run(
        """
        image = readfits('{{ path }}/image.fits')
        psf = readfits('{{ path }}/psf.fits')
        out = aia_deconvolve_richardsonlucy(image, psf, niter=25)
        """,
        args={"path": tmp_path.as_posix()},
        save_vars=["out"],
        verbose=False,
    )
    map_decon = aiapy.psf.deconvolve(aia_171_positive_map, psf=psf_128, iterations=25)
    np.testing.assert_allclose(map_decon.data, ssw["out"], rtol=1e-10)


@pytest.mark.remote_data
def test_register_consistent(idl_environment, aia_193_level1_map, tmp_path) -> None:
    """
    `register` puts a full-resolution image on the same grid as SSW ``aia_prep``.

    SSW interpolates with damped cubic convolution and stores integer DN, so pixels
    differ by about 0.5% on disk; a 0.25 pixel misregistration doubles that.
    """
    filename = tmp_path / "aia_lev1.fits"
    with pytest.warns(VerifyWarning, match="The 'BLANK' keyword is only applicable to integer data"):
        aia_193_level1_map.save(filename)
    ssw = idl_environment.run(
        """
        mreadfits, '{{ filename }}', index, data
        aia_prep, index, data, oindex, odata, /no_mpo_update
        crpix = [oindex.crpix1, oindex.crpix2]
        cdelt = [oindex.cdelt1, oindex.cdelt2]
        """,
        save_vars=["odata", "crpix", "cdelt"],
        args={"filename": filename.as_posix()},
        verbose=False,
    )
    registered = register(aia_193_level1_map)
    np.testing.assert_allclose([registered.meta["crpix1"], registered.meta["crpix2"]], ssw["crpix"])
    assert_quantity_allclose(registered.scale, ssw["cdelt"] * u.arcsec / u.pix)
    assert registered.data.shape == ssw["odata"].shape
    y, x = np.indices(registered.data.shape)
    disk = (
        np.hypot(x - registered.meta["crpix1"] + 1, y - registered.meta["crpix2"] + 1) < 0.95 * registered.meta["r_sun"]
    )
    difference = np.abs(registered.data - ssw["odata"])[disk]
    assert np.median(difference) / np.median(ssw["odata"][disk]) < 0.0075
