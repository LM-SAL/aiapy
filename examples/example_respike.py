"""
========================================
Re-spiking Level 1 Data
========================================

This example demonstrates how to "re-spike" AIA level 1 images. AIA level 1
images have been corrected for hot-pixels (commonly referred to as "spikes")
using an automated correction algorith, which detects them, removes them and
replaces the "holes" left in the image via interpolation. However, For certain
research topics, this automated hot-pixel removal process may result in unwanted
removal of bright points which may be real. In this example, we will demonstrate
how to do revert this removal by putting back all the removed pixel values with
`respike()` in `aiapy`.
This corresponds to the `aia_respike.pro` IDL procedure as described in the
`SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`_.

"""

import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from sunpy.net import Fido, attrs
from astropy.io import fits
import drms
import sunpy.map
from aiapy.calibrate import respike

###########################################################
#The keywords (`LVL_NUM`, `NSPIKES`) describe the level number of the AIA data
#(i.e. level 1) and how many hot pixels were removed from the image (i.e.
#the "spikes").
#The data containing the information of the pixel position and the intensities
#of the hot pixels  are available as a separate segment of the
#`aia.lev1_euv_12s` and `aia.lev1_uv_24s` data series from the
#`Joint Science Operations Center (JSOC) <http://jsoc.stanford.edu/>`_.
#
#Below we show how this data can be retrieved.
#
#First, let's fetch a level 1 AIA image and its associated spike data from the `JSOC <https://jsoc.stanford.edu/>`_. For our demonstration we use a 193 Å image taken on 15 March 2013 and create a `~sunpy.map.Map` object.

q = Fido.search(attrs.Time('2013-03-15T12:01:00', '2013-03-15T12:01:10'),
    attrs.jsoc.Wavelength(193*u.angstrom),
    attrs.jsoc.Series('aia.lev1_euv_12s'),
    attrs.jsoc.Segment('image') & attrs.jsoc.Segment('spikes'),
    attrs.jsoc.Notify('sunpy@sunpy.org'))

ms = Fido.fetch(q)

#ensure that the query returns 'image' filename as the first entry in the list and 'spikes' as the second
sortms = list(ms)
for i in range(len(sortms)): sortms[i] = [s for s in ms if ['image','spikes'][i] in s]
ms = [i[0] for i in sortms]

#read fetched image and produce a sunpy map from the downloaded level 1 193 A snapshot
dummy, image = fits.open(ms[0])
image.verify('fix')
AIAMap1 = sunpy.map.Map(image.data, image.header)

#read the data array for the spikes (note that there is no useful header information for spikes)
dummy, spikes = fits.open(ms[1])
spikes.verify('fix')
spikes = spikes.data


###########################################################
#The spike data are stored as separate data segments in JSOC
#as a [3, N] arrays, where N the number of spikes removed.
#
#[0,:] contains the 1-D pixel index representing its (y, x) position
#in the 2-D image data,
#[1,:] contains the intensity value for each spike removed, and
#[2,:] contains the intensity value used in replacing the spike removed (via interpolation).
#
#Spike data come without header information - the pixel positions are
#given with respect to the level 0.5 full disk image.
#
#The `respike()` function fetches the spike data for a given AIAmap from JSOC
#However, it is advised that the user should get the data (e.g., as shown above)
#and pass the spike array to `respike()` via the arguments, i.e.,
#`respike(AIAMap, spikes=spike_data_array)`
#

###########################################################
#Now we are ready to respike the Level 1 AIA image. We produce a sunpy
#AIAMap and run `respike(AIAMap, spikes=spikes)`.
#The routine spike_coords() also exports the pixel x and y positions as
#
#`spike_pixcoords = spike_coords(AIAMap, spikes=spikes, pixels=True)`

#perform respike operation on the desired map.
#Note that this operation alters both the data and the metadata (the latter by updating the
#keywords LVL_NUM, NSPIKES and COMMENTS). Using `spike_coords()` one can get the pixel position of the spikes
#in units (arcsec) or in pixels (if pixels=True).
AIAMap = respike(AIAMap1, spikes = spikes)

#export positions of spikes in units (arcsec)
hpc_spikes = spike_coords(AIAMap1, spikes = spikes)

###########################################################
#Alternatively, given an AIAMap, `respike()` can fetch the associated spike data from JSOC.
#This can be done by simply typing `AIAMap=respike(AIAMap1)`
#
#Tip: Fetching the spike data automatically through the functions seems reasonable
#if the number of AIA maps in the image series is small. However,
#if the number of AIA maps is large, it is more efficient to request the spike data
#for the image series via JSOC and then import them using
#the `spike=` argument in these functions
#

#Now let's plot a cutout of the 'clean' AIAMap1 (Level 1) and "re-spiked" AIAMap (i.e. Level 0.5) for a region with hot pixels
top_right = SkyCoord(30 * u.arcsec, 420 * u.arcsec, frame = AIAMap1.coordinate_frame)
bottom_left = SkyCoord(-120 * u.arcsec, 280 * u.arcsec, frame = AIAMap1.coordinate_frame)
submap = AIAMap.submap(bottom_left, top_right)
submap1 = AIAMap1.submap(bottom_left, top_right)

#plot the two maps for comparison
#Level 1
fig = plt.figure()
ax = fig.add_subplot(111, projection = submap)
ax.plot_coord(hpc_spikes, 'o', color = 'white', fillstyle = 'None', markersize = 15)
fig.subplots_adjust(top = 0.87)
fig.suptitle('Level 1 data ("De-spiked")')
submap1.plot()
plt.show()

#Level 0.5
fig = plt.figure()
ax = fig.add_subplot(111, projection = submap1)
ax.plot_coord(hpc_spikes, 'o', color = 'white', fillstyle = 'None', markersize = 15)
fig.subplots_adjust(top = 0.87)
fig.suptitle('"Re-spiked" data (back to Level 0.5)')
submap.plot()
plt.show()

###########################################################
#The white circles show the location of the spikes. We plot the circles in both
#Level 1 and Level 0.5 data to aid the visual inspection.
#
#Lastly, let's check the metadata information:

print("Data downloaded from VSO/JSOC are Level {} and had {} hot pixels removed".format(submap1.meta['lvl_num'], submap1.meta['nspikes']))
print("These data when respiked are converted to Level {} with {} hot pixels removed (i.e. all hotpixels are now included)".format(submap.meta['lvl_num'], submap.meta['nspikes']))

###########################################################
#This confirms that the data have been successfuly converted to Level 0.5
#after re-spiking and the header information is also reflecting this
#conversion. You may now proceed and use custom methods for the removal of
#hot pixels from such Level 0.5 AIA maps.
