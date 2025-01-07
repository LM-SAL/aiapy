"""
=============================================
Updating pointing keywords in the FITS header
=============================================

This example demonstrates how to update the metadata in
an AIA FITS file to ensure that it has the most accurate
information regarding the spacecraft pointing.
"""

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map

import aiapy.data.sample as sample_data
from aiapy.calibrate import update_pointing
from aiapy.calibrate.util import get_pointing_table

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
# `aiapy.calibrate.update_pointing` function.
# One needs to get the pointing information from the JSOC using the
# `aiapy.calibrate.util.get_pointing_table` function first.

# Make range wide enough to get closest 3-hour pointing
pointing_table = get_pointing_table("JSOC", start=aia_map.date - 12 * u.h, end=aia_map.date + 12 * u.h)
aia_map_updated_pointing = update_pointing(aia_map, pointing_table=pointing_table)

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
# location in the header are correct.

print(aia_map_updated_pointing.meta["hgln_obs"])
print(aia_map_updated_pointing.meta["hglt_obs"])

###############################################################################
# Finally, plot the fixed map.

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map_updated_pointing)
aia_map_updated_pointing.plot(axes=ax)

plt.show()
