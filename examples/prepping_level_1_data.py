"""
=====================================
Registering and aligning level 1 data
=====================================

This example demonstrates how to convert AIA images to a common pointing,
rescale them to a common plate scale, and remove the roll angle.
This process is often referred to as "aia_prep" and the resulting data are typically referred to as level 1.5 data.
In this example, we will demonstrate how to do this with `aiapy`.
This corresponds to the ``aia_prep.pro`` procedure as described in the `SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`__.
"""

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map

import aiapy.data.sample as sample_data
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.util import get_pointing_table

###############################################################################
# Performing multi-wavelength analysis on level 1 data can be problematic as
# each of the AIA channels have slightly different spatial scales and roll
# angles. Furthermore, the estimates of the pointing keywords (``CDELT1``, ``CDELT2``, ``CRPIX1``,
# ``CRPIX2``, ``CROTA2``) may have been improved due to limb fitting procedures after
# the level 1 file has been created.
# The `Joint Science Operations Center (JSOC) <http://jsoc.stanford.edu/>`_ stores
# AIA image data and metadata separately; when users download AIA data, these
# two data types are combined to produce a FITS file. While metadata are
# continuously updated at the JSOC, previously downloaded FITS files will not
# contain the most recent information.
#
# Thus, before performing any multi-wavelength analyses, level 1 data
# should be updated to the most recent and accurate pointing and interpolated
# to a common grid in which the y-axis of the image is aligned
# with solar North.
#
# First, let's read a level 1 94 Å AIA image from the ``aiapy`` sample data into
# a `~sunpy.map.Map` object.

aia_map = sunpy.map.Map(sample_data.AIA_094_IMAGE)

###############################################################################
# The first step in this process is to update the metadata of the map to the
# most recent pointing using  the `aiapy.calibrate.update_pointing` function.
# One needs to get the pointing information from the JSOC using the
# `aiapy.calibrate.util.get_pointing_table` function.

# Make range wide enough to get closest 3-hour pointing
pointing_table = get_pointing_table("JSOC", time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h))
aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)

###############################################################################
# If we take a look at the plate scale and rotation matrix of the map, we
# find that the scale is slightly off from the expected value of :math:`0.6''` per
# pixel and that the rotation matrix has off-diagonal entries.

print(aia_map_updated_pointing.scale)
print(aia_map_updated_pointing.rotation_matrix)

###############################################################################
# We can use the `aiapy.calibrate.register` function to scale the image to
# the :math:`0.6''` per pixel and derotate the image such that the y-axis is aligned
# with solar North.

aia_map_registered = register(aia_map_updated_pointing)

###############################################################################
# If we look again at the plate scale and rotation matrix, we
# should find that the plate scale in each direction is :math:`0.6''`
# per pixel and that the rotation matrix is diagonalized.
# The image in ``aia_map_registered`` is now a level 1.5 data product.

print(aia_map_registered.scale)
print(aia_map_registered.rotation_matrix)

###############################################################################
# Finally, we can plot the exposure-normalized map.
# Note that small negative pixel values are possible because
# CCD images were taken with a pedestal set at ~100 DN.
# This pedestal is then subtracted when the JSOC pipeline
# performs dark (+ pedestal) subtraction and flat fielding
# to generate level 1 files.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map_registered)
aia_map_registered.plot(axes=ax)

plt.show()
