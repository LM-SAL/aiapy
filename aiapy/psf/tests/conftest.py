"""
Shared fixtures for PSF tests.
"""

import pytest

import aiapy.psf


@pytest.fixture()
def psf(channels):
    return aiapy.psf.psf(channels[0], use_preflightcore=True, diffraction_orders=[-1, 0, 1])
