"""
Shared fixtures for PSF tests.
"""

import pytest

import aiapy.psf
from aiapy.conftest import CHANNELS


@pytest.fixture
def psf():
    return aiapy.psf.psf(CHANNELS[0], use_preflightcore=True, diffraction_orders=[-1, 0, 1])
