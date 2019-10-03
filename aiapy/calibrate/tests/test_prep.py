"""
Tests for functions that calibrate/prep AIA image data
"""
import tempfile

import numpy as np
import pytest
from astropy.io.fits.verify import VerifyWarning
import sunpy.data.test
from sunpy.map import Map

from aiapy.calibrate import register


@pytest.fixture
def original():
    return Map(sunpy.data.test.get_test_filepath('aia_171_level1.fits'))


@pytest.fixture
def lvl_15_map(original):
    return register(original)


def test_register(original, lvl_15_map):
    """
    Test that header info for the map has been correctly updated after the
    map has been scaled to 0.6 arcsec / pixel and aligned with solar north
    """
    # Check all of these for Map attributes and .meta values?
    # Check array shape
    assert lvl_15_map.data.shape == original.data.shape
    # Check crpix values
    assert lvl_15_map.meta['crpix1'] == lvl_15_map.data.shape[1] / 2.0 + 0.5
    assert lvl_15_map.meta['crpix2'] == lvl_15_map.data.shape[0] / 2.0 + 0.5
    # Check cdelt values
    assert lvl_15_map.meta['cdelt1'] / 0.6 == int(lvl_15_map.meta['cdelt1'] / 0.6)
    assert lvl_15_map.meta['cdelt2'] / 0.6 == int(lvl_15_map.meta['cdelt2'] / 0.6)
    # Check rotation value, I am assuming that the inaccuracy in
    # the CROTA -> PCi_j matrix is causing the inaccuracy here
    np.testing.assert_allclose(
        lvl_15_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert lvl_15_map.meta['lvl_num'] == 1.5


def test_filesave(lvl_15_map):
    """
    Test that adjusted header values are still correct after saving the map
    and reloading it.
    """
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    with pytest.warns(
            VerifyWarning,
            match="The 'BLANK' keyword is only applicable to integer data"):
        lvl_15_map.save(afilename, overwrite=True)
    load_map = Map(afilename)
    # Check crpix values
    assert load_map.meta['crpix1'] == lvl_15_map.data.shape[1] / 2.0 + 0.5
    assert load_map.meta['crpix2'] == lvl_15_map.data.shape[0] / 2.0 + 0.5
    # Check cdelt values
    assert load_map.meta['cdelt1'] / 0.6 == int(load_map.meta['cdelt1'] / 0.6)
    assert load_map.meta['cdelt2'] / 0.6 == int(load_map.meta['cdelt2'] / 0.6)
    # Check rotation value
    np.testing.assert_allclose(
        lvl_15_map.rotation_matrix, np.identity(2), rtol=1e-5, atol=1e-8)
    # Check level number
    assert load_map.meta['lvl_num'] == 1.5


def test_unsupported_maps(original):
    """
    Make sure we raise an error when an unsupported map is passed in
    """
    # A submap
    original_cutout = original.submap(original.center,
                                      original.top_right_coord)
    with pytest.raises(ValueError):
        _ = register(original_cutout)
    # A Map besides AIA or HMI
    non_sdo_map = Map(sunpy.data.test.get_test_filepath(
        'mdi_fd_Ic_6h_01d.5871.0000_s.fits'))
    with pytest.raises(ValueError):
        _ = register(non_sdo_map)
