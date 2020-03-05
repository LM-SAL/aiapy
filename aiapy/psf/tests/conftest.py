"""
Shared fixtures for PSF tests
"""
import pytest
import astropy.units as u

import aiapy.psf


@pytest.fixture(scope='module')
def channels():
    return [94, 131, 171, 193, 211, 304, 335] * u.angstrom


@pytest.fixture
def psf(channels):
    return aiapy.psf.psf(channels[0],
                         use_preflightcore=True,
                         diffraction_orders=[-1, 0, 1])
