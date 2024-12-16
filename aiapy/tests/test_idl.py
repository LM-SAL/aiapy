"""
Contains all the the IDL specific tests for aiapy.
"""

from pathlib import Path

import numpy as np
import pytest

import astropy.units as u

from aiapy.calibrate import estimate_error
from aiapy.conftest import CHANNELS


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
    error_table = Path(idl_environment.ssw_home) / "sdo" / "aia" / "response" / "aia_V3_error_table.txt"
    ssw = idl_environment.run(
        idl,
        save_vars=["error"],
        args={
            "channel": channel.to("angstrom").value,
            "data": counts.to("ct pixel-1").value,
            "error_table": error_table,
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
    assert u.allclose(error, error_ssw, rtol=1e-4)


@pytest.fixture(params=CHANNELS)
def psf_idl(idl_environment, request):
    """
    The point spread function as calculated by aia_calc_psf.pro.
    """
    r = idl_environment.run(
        "psf = aia_calc_psf({{channel}},/use_preflightcore)",
        args={"channel": f"{request.value:.0f}"},
        save_vars=["psf"],
        verbose=False,
    )
    return r["psf"]


def test_psf_consistent(psf_94, psf_idl) -> None:
    """
    Check whether PSF is consistent with IDL calculation.

    **WARNING** This test will take a very long time to run.
    """
    # NOTE: The IDL and Python PSF functions have been found to
    # agree within 0.2% for all points along the PSF arms for
    # both the preflight and non-preflight cases.
    # NOTE: Only compare values above some threshold as the
    # rest of the PSF is essentially noise
    i_valid = np.where(psf_idl > 1e-10)
    assert np.allclose(psf_94[i_valid], psf_idl[i_valid], atol=0.0, rtol=2e-3)
