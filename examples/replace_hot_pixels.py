"""
========================================
Re-spiking Level 1 Images
========================================

This example demonstrates how to "re-spike" AIA level 1 images
"""

import matplotlib.pyplot as plt

import astropy.units as u
import sunpy.map
from astropy.coordinates import SkyCoord
from sunpy.net import Fido, attrs

from aiapy.calibrate import fetch_spikes, respike

####################################################
# AIA level 1 images have been corrected for hot-pixels (commonly referred to
# as "spikes") using an automated correction algorithm which detects them,
# removes them, and replaces the "holes" left in the image via interpolation.
# However, for certain research topics, this automated hot-pixel removal
# process may result in unwanted removal of bright points which may be
# physically meaningful. In this example, we will demonstrate how to revert
# this removal by putting back all the removed pixel values with the
# `~aiapy.calibrate.respike` in function. This corresponds to the
# `aia_respike.pro` IDL procedure as described in the
# `SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`_.
#
# The header keywords `LVL_NUM` and `NSPIKES` describe the level number of the
# AIA data (e.g. level 1) and how many hot pixels were removed from the image
# (i.e. the "spikes"). The data containing the information of the pixel
# position and the intensities of the removed hot pixels are available from the
# `Joint Science Operations Center (JSOC) <http://jsoc.stanford.edu/>`_ as a
# separate segment of the `aia.lev1_euv_12s` and `aia.lev1_uv_24s` data series

####################################################
# First, let's fetch a level 1 AIA image and read it into a `~sunpy.map.Map`. For
# our demonstration, we use a 193 Å image taken on 15 March 2013.
q = Fido.search(attrs.Time('2013-03-15T12:01:00', '2013-03-15T12:01:10'),
                attrs.Wavelength(193*u.angstrom),
                attrs.Instrument('AIA'))
f = Fido.fetch(q)
m = sunpy.map.Map(f)


###########################################################
# The spike data are stored as separate data segments in JSOC
# as a :math:`3\times N` arrays, where :math:`N` is the number of spikes
# removed and the three dimensions correspond to the the 1-D pixel index
# of the spike, intensity value of the removed spikes, and the intensity value
# used in replacing the removed spike (via interpolation).
# The spike pixel positions are given with respect to the level 1 full-disk
# image.
#
# We can use the `~aiapy.calibrate.fetch_spikes` function to query the JSOC
# for the spike positions and intensity values and convert the positions of the
# spikes to the 2D pixel full-disk pixel coordinate system given a
# `~sunpy.map.Map` representing a level 1 AIA image.
#
positions, values = fetch_spikes(m)

###########################################################
# Now we are ready to respike the level 1 AIA image. The
# `~aiapy.calibrate.respike` function performs the respike operation on the given
# input image and returns a `~sunpy.map.Map` with the respiked image. This
# operation also alters the metadata by updating the `LVL_NUM`, `NSPIKES`,
# and `COMMENTS` keywords.
#
# Note that explicitly specifying the spike positions and values is optional.
# If they are not given, they are automatically queried from the JSOC.
m_respiked = respike(m, spikes=(positions, values))

###########################################################
# Now let's create a cutouts of the original level 1 and "re-spiked" (i.e.
# level 0.5) images for a region with hot pixels.
top_right = SkyCoord(30 * u.arcsec, 420 * u.arcsec,
                     frame=m.coordinate_frame)
bottom_left = SkyCoord(-120 * u.arcsec, 280 * u.arcsec,
                       frame=m.coordinate_frame)
m_cutout = m.submap(bottom_left, top_right=top_right)
m_respiked_cutout = m_respiked.submap(bottom_left, top_right=top_right)

###########################################################
# Note that we can also retrieve the positions of the spikes
# as `~astropy.coordinates.SkyCoord` objects in the projected coordinate
# system of the image using the `as_coords=True` keyword argument. This
# gives us only those spikes in the field of view of the cutout.
spike_coords, _ = fetch_spikes(m_cutout, as_coords=True)

###########################################################
# Finally, let's plot the two cutouts for comparison and plot
# the positions of the spikes in both images, denoted by white
# circles.
fig = plt.figure()
ax = fig.add_subplot(121, projection=m_cutout)
ax.plot_coord(spike_coords, 'o', color='white', fillstyle='none',
              markersize=15)
m_cutout.plot(axes=ax, title='Level 1 "de-spiked" data')
lon, lat = ax.coords
lon.set_axislabel('HPC Longitude')
lat.set_axislabel('HPC Latitude')
ax = fig.add_subplot(122, projection=m_respiked_cutout)
ax.plot_coord(spike_coords, 'o', color='white', fillstyle='none',
              markersize=15)
m_respiked_cutout.plot(axes=ax, annotate=False)
ax.set_title('Level 0.5 "re-spiked" data')
lon, lat = ax.coords
lon.set_axislabel('HPC Longitude')
lat.set_axislabel(' ')
lat.set_ticklabel_visible(False)
plt.show()

###########################################################
# Lastly, let's check the metadata in both the level 1 and resulting
# 0.5 images to double check that the appropriate keywords have been updated.
for k in ['lvl_num', 'nspikes', 'comments']:
    print(f'Level 1: {k}: {m_cutout.meta.get(k)}')
    print(f'Level 0.5: {k}: {m_respiked_cutout.meta.get(k)}')
