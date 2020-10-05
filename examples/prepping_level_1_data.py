"""
=======================================
Registering and Aligning Level 1 Data
=======================================

This example demonstrates how to convert AIA images to a common pointing,
rescale them to a common plate scale, and remove the roll angle. This process
is often referred to as "aia_prep" and the resulting data are typically
referred to as level 1.5 data. In this example, we will demonstrate how to do
this with `aiapy`. This corresponds to the `aia_prep.pro` procedure as
described in the `SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`_.
"""

import astropy.units as u
from sunpy.net import Fido, attrs
import sunpy.map

from aiapy.calibrate import register, update_pointing, normalize_exposure

###########################################################
# Performing multi-wavelength analysis on level 1 data can be problematic as
# each of the AIA channels have slightly different spatial scales and roll
# angles. Furthermore, the estimates of the pointing keywords (`CDELT1`, `CDELT2`, `CRPIX1`,
# `CRPIX2`, `CROTA2`) may have been improved due to limb fitting procedures. The
# `Joint Science Operations Center (JSOC) <http://jsoc.stanford.edu/>`_ stores
# AIA image data and metadata separately; when users download AIA data, these
# two data types are combined to produce a FITS file. While metadata are
# continuously updated at JSOC, previously downloaded FITS files will not
# contain the most recent information.
#
# Thus, before performing any multi-wavelength analyses, level 1 data
# should be updated to the most recent and accurate pointing and interpolated
# to a common grid in which the y-axis of the image is aligned
# with solar North.
#
# First, let's fetch level 1 AIA images from the
# `Virtual Solar Observatory <https://sdac.virtualsolar.org/cgi/search>`_
# from 1 January 2019 for the 94 Å channel and create a `~sunpy.map.Map`
# object.
q = Fido.search(
    attrs.Time('2019-01-01T00:00:00', '2019-01-01T00:00:11'),
    attrs.Instrument('AIA'),
    attrs.Wavelength(wavemin=94*u.angstrom, wavemax=94*u.angstrom),
)
m = sunpy.map.Map(Fido.fetch(q))

###########################################################
# The first step in this process is to update the metadata of the map to the
# most recent pointing using  the `~aiapy.calibrate.update_pointing` function.
# This function queries the JSOC for the most recent pointing information,
# updates the metadata, and returns a `~sunpy.map.Map` with updated metadata.
m_updated_pointing = update_pointing(m)

###########################################################
# If we take a look at the plate scale and rotation matrix of the map, we
# find that the scale is slightly off from the expected value of :math:`0.6''` per
# pixel and that the rotation matrix has off-diagonal entries.
print(m_updated_pointing.scale)
print(m_updated_pointing.rotation_matrix)

###########################################################
# We can use the `~aiapy.calibrate.register` function to scale the image to
# the :math:`0.6''` per pixel and derotate the image such that the y-axis is aligned
# with solar North.
m_registered = register(m_updated_pointing)

###########################################################
# If we look again at the plate scale and rotation matrix, we
# should find that the plate scale in each direction is :math:`0.6''`
# per pixel and that the rotation matrix is diagonalized.
# The image in `m_registered` is now a level 1.5 data product.
print(m_registered.scale)
print(m_registered.rotation_matrix)

###########################################################
# Though it is not typically part of the level 1.5 "prep" data pipeline,
# it is also common to normalize the image to the exposure time such that
# the units of the image are DN / pixel / s.
m_normalized = normalize_exposure(m_registered)

###########################################################
# Finally, we can plot the exposure-normalized map.
# Note that small negative pixel values are possible because
# CCD images were taken with a pedestal set at ~ 100 DN.
# This pedestal is then subtracted when the JSOC pipeline
# performs dark (+pedestal) subtraction and flatfielding
# to generate level 1 files.
m_normalized.peek(vmin=0)
