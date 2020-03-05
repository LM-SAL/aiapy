"""
Test PSF calculation
"""
import numpy as np
import pytest

import aiapy.psf

MESH_PROPERTIES = [
    'angle_arm',
    'error_angle_arm',
    'spacing_e',
    'spacing_fp',
    'mesh_pitch',
    'mesh_width',
    'width',
]


@pytest.fixture
def psf_full(channels):
    return aiapy.psf.psf(channels[0], use_preflightcore=True)


@pytest.fixture(scope='module')
def psf_idl(idl_environment, channels):
    """
    The point spread function as calculated by `aia_calc_psf.pro`
    """
    r = idl_environment.run('psf = aia_calc_psf({{channel}},/use_preflightcore)',
                            args={'channel': f'{channels[0].value:.0f}'},
                            save_vars=['psf'],
                            verbose=False)
    return r['psf']


@pytest.mark.parametrize('use_preflightcore', [True, False])
def test_filter_mesh_parameters(use_preflightcore, channels):
    params = aiapy.psf.filter_mesh_parameters(
        use_preflightcore=use_preflightcore)
    assert isinstance(params, dict)
    assert all([c in params for c in channels])
    assert all([all([p in params[c] for p in MESH_PROPERTIES])
                for c in channels])


def test_psf(psf):
    assert isinstance(psf, np.ndarray)
    assert psf.shape == (4096, 4096)


def test_psf_consistent(psf_full, psf_idl):
    """
    Check whether PSF is consistent with IDL calculation.

    .. note:: This test will take a very long time to run.
    """
    # NOTE: The IDL and Python PSF functions have been found to
    # agree within 0.2% for all points along the PSF arms for
    # both the preflight and non-preflight cases.
    # NOTE: Only compare values above some threshold as the
    # rest of the PSF is essentially noise
    i_valid = np.where(psf_idl > 1e-10)
    assert np.allclose(psf_full[i_valid], psf_idl[i_valid], atol=0., rtol=2e-3)
