"""
==========================================================
Updating pointing and observer keywords in the FITS header
==========================================================

This example demonstrates how to update the metadata in
an AIA FITS file to ensure that it has the most accurate
information regarding the spacecraft pointing and observer
position.
"""

import matplotlib.pyplot as plt
import sunpy.map

import aiapy.data.sample as sample_data
from aiapy.calibrate import fix_observer_location, update_pointing

###############################################################################
# An AIA FITS header contains various pieces of
# `standard <https://fits.gsfc.nasa.gov/fits_standard.html>`_.
# metadata that are critical to the physical interpretation of the data.
# These include the pointing of the spacecraft, necessary for connecting
# positions on the pixel grid to physical locations on the Sun, as well as
# the observer (i.e., satellite) location.
#
# While this metadata is recorded in the FITS header, some values in
# the headers exported by data providers (e.g.
# `Joint Science Operations Center (JSOC) <http://jsoc.stanford.edu/>`_ and
# the `Virtual Solar Observatory <https://sdac.virtualsolar.org/cgi/search>`_
# may not always be the most accurate. In the case of the spacecraft
# pointing, a more accurate 3-hourly pointing table is available from the
# JSOC.
#
# For this example, we will read a 171 Å image from the aiapy sample data
# into a `~sunpy.map.Map` object.

aia_map = sunpy.map.Map(sample_data.AIA_171_IMAGE)

###############################################################################
# To update the pointing keywords, we can pass our `~sunpy.map.Map` to the
# `aiapy.calibrate.update_pointing` function. This function will query the
# JSOC, using `~sunpy`, for the most recent pointing information, update
# the metadata, and then return a new `~sunpy.map.Map` with this updated
# metadata.

aia_map_updated_pointing = update_pointing(aia_map)

###############################################################################
# If we inspect the reference pixel and rotation matrix of the original map:

print(aia_map.reference_pixel)
print(aia_map.rotation_matrix)

###############################################################################
# and the map with the updated pointing information:

print(aia_map_updated_pointing.reference_pixel)
print(aia_map_updated_pointing.rotation_matrix)

###############################################################################
# We find that the relevant keywords, ``CRPIX1``, ``CRPIX2``, ``CDELT1``, ``CDELT2``,
# and ``CROTA2``, have been updated.
#
# Similarly, the Heliographic Stonyhurst (HGS) coordinates of the observer
# location in the header are inaccurate. If we check the HGS longitude keyword
# in the header, we find that it is 0 degrees which is not the HGS longitude
# coordinate of SDO.

print(aia_map_updated_pointing.meta["hgln_obs"])
print(aia_map_updated_pointing.meta["hglt_obs"])

###############################################################################
# To update the HGS observer coordinates, we can use the
# `aiapy.calibrate.fix_observer_location` function. This function reads the
# correct observer location from Heliocentric Aries Ecliptic (HAE) coordinates
# in the header, converts them to HGS, and replaces the inaccurate HGS
# keywords.

aia_map_observer_fixed = fix_observer_location(aia_map_updated_pointing)

###############################################################################
# Looking again at the HGS longitude and latitude keywords, we can see that
# they have been updated.
print(aia_map_observer_fixed.meta["hgln_obs"])
print(aia_map_observer_fixed.meta["hglt_obs"])

###############################################################################
# Note that in `~sunpy.map.sources.AIAMap`, the `~sunpy.map.GenericMap.observer_coordinate`
# attribute is already derived from the HAE coordinates such that it is not
# strictly necessary to apply `aiapy.calibrate.fix_observer_location`. For
# example, the unfixed `~sunpy.map.Map` will still have an accurate derived
# observer position

print(aia_map_updated_pointing.observer_coordinate)

###############################################################################
# However, we suggest that users apply this fix such that the information
# stored in `~sunpy.map.GenericMap.meta` is accurate and consistent.
#
# Finally, plot the fixed map.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map_observer_fixed)
aia_map_observer_fixed.plot(axes=ax)

plt.show()
