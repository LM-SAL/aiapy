"""
============================================
Requesting specific AIA images from the JSOC
============================================

This example shows how to request a specific series of AIA images from the JSOC.

We will be filtering the data we require by keywords and requesting short exposure images from a recent (at time of writing) flare.
"""

import os

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import AsinhStretch, ImageNormalize

import sunpy.map
from sunpy.net import Fido, attrs

from aiapy.calibrate import correct_degradation, register, update_pointing
from aiapy.calibrate.util import get_correction_table, get_pointing_table

###############################################################################
# Exporting data from the JSOC requires registering your
# email first. Please replace the text after the ``=``
# with your email address once you have registered.
# `See this page for more details. <http://jsoc.stanford.edu/ajax/register_email.html>`__

jsoc_email = os.environ.get("JSOC_EMAIL")

###############################################################################
# Our goal is to request data of a recent X-class flare.
# The X-class flare occurred on the 2021/07/03 at 14:30:00 UTC.
# We will focus on the 5 minutes before and after this time.
# What we want to do is only get the shorter exposures,
# which are going to be the flare related.

query = Fido.search(
    attrs.Time("2021-07-03 14:25:00", "2021-07-03 14:35:00"),
    attrs.jsoc.Series("aia.lev1_euv_12s"),
    attrs.Wavelength(211 * u.AA),
    attrs.jsoc.Notify(jsoc_email),
    attrs.jsoc.Keyword("EXPTIME") <= 2,
    attrs.jsoc.Segment("image"),
)

print(query)

###############################################################################
# Now we will download the data and "aia prep" the
# data with every feature of `aiapy` and plot the
# data sequence using `sunpy`.

files = Fido.fetch(query)

level_1_maps = sunpy.map.Map(files)
# We get the pointing table outside of the loop for the relevant time range.
# Otherwise you're making a call to the JSOC every single time.
pointing_table = get_pointing_table(
    "jsoc", time_range=(level_1_maps[0].date - 3 * u.h, level_1_maps[-1].date + 3 * u.h)
)
# The same applies for the correction table.
correction_table = get_correction_table(source="jsoc")

level_15_maps = []
for a_map in level_1_maps:
    map_updated_pointing = update_pointing(a_map, pointing_table=pointing_table)
    map_registered = register(map_updated_pointing)
    map_degradation = correct_degradation(map_registered, correction_table=correction_table)
    map_normalized = map_degradation / map_degradation.exposure_time
    bottom_left = SkyCoord(500 * u.arcsec, 100 * u.arcsec, frame=map_normalized.coordinate_frame)
    top_right = SkyCoord(1500 * u.arcsec, 700 * u.arcsec, frame=map_normalized.coordinate_frame)
    map_cropped = map_normalized.submap(bottom_left, top_right=top_right)
    level_15_maps.append(map_cropped)

###############################################################################
# Finally, we create a sequence of maps and animate it:

sequence = sunpy.map.Map(level_15_maps, sequence=True)

fig = plt.figure()
ax = fig.add_subplot(projection=sequence.maps[0])
ani = sequence.plot(axes=ax, norm=ImageNormalize(vmin=0, vmax=1e3, stretch=AsinhStretch()))

plt.show()
