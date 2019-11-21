"""
=======================================
Registering and Aligning Level 1 Data
=======================================

This example demonstrates how to convert AIA images
to a common pointing, rescale them to a common resolution,
and remove the roll angle. This corresponds to the `aia_prep.pro`
procedure as described in the `SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`_
"""

import astropy.units as u
from sunpy.net import Fido, attrs
import sunpy.map

from aiapy.calibrate import register, update_pointing

###########################################################
# Performing multi-wavelength analysis on level 1 data can be problematic as
# each of the AIA channels have slightly different spatial scales and roll
# angles. Furthermore, the pointing keywords (`CDELT1`, `CDELT2`, `CRPIX1`,
# `CRPIX2`, `CROTA2`) may be out of date. Thus, before performing any
# multi-wavelength analyses, level 1 data should be updated to the most
# recent and accurate pointing and interpolated to a common resolution
# grid in which the y-axis of the image is aligned with solar North.
#
# This process is often referred to as "aia_prep" and the resulting data are
# typically referred to as level 1.5 data. In this example, we will
# demonstrate how to do this with `aiapy`.
#
# First, let's fetch level 1 AIA images from the VSO from 1 January
# 2019 for the 94 Å and 131 Å channels and create a `~sunpy.map.Map` object.
q = Fido.search(
    attrs.Time('2019-01-01T00:00:00', '2019-01-01T00:00:11'),
    attrs.Instrument('AIA'),
    attrs.Wavelength(wavemin=94*u.angstrom, wavemax=131*u.angstrom),
)
maps = sunpy.map.Map(Fido.fetch(q))

###########################################################
# The first step in this process is to update the metadata of the map to the
# most recent pointing using  the `~aiapy.calibrate.update_pointing` function.
# This function queries the JSOC for the most recent pointing information,
# updates the metadata, and returns a `~sunpy.map.Map` with updated metadata.
maps_updated_pointing = [update_pointing(m) for m in maps]

###########################################################
# If we take a look at the plate scale and rotation matrices of each map, we
# find that each scale is slightly different for each channel and that the
# rotation matrices have off-diagonal entries.
for m in maps_updated_pointing:
    print(m.scale)
for m in maps_updated_pointing:
    print(m.rotation_matrix)

###########################################################
# We can use the `~aiapy.calibrate.register` function to scale each image to
# a common resolution and derotate the image such that the y-axis of the image
# is aligned with solar North. To save memory and compute time, we will first
# resample each map to have 1024 pixels in each direction. This step can be
# omitted if you'd prefer to work with the full-resolution images.
maps_resampled = [m.resample((1024, 1024)*u.pixel) for m in maps_updated_pointing]
maps_registered = [register(m) for m in maps_updated_pointing]

###########################################################
# If we look again at the plate scale and rotation matrix of each image, we
# should find that the resolution in each direction of all images is
# the same and that the rotation matrices are all diagonalized.
# The images in `maps_registered` are now level 1.5 data products. Note that
# if we had not resampled the image prior to registration, the new plate scale
# of each image would be 0.6 arcseconds per pixel.
for m in maps_registered:
    print(m.scale)
for m in maps_registered:
    print(m.rotation_matrix)

###########################################################
# Though it is not typically part of the level 1.5 "prep" data pipeline,
# it is also common to normalize the image to the exposure time such that
# the units of the image are DN / pixel / s. This is straightforward to do
# using the information exposed by the `~sunpy.map.Map` API.
maps_normalized = [sunpy.map.Map(m.data/m.exposure_time.to(u.s).value, m.meta)
                   for m in maps_registered]
