"""
===================================================
Isolating Fe XVIII emission in the AIA 94 Å channel
===================================================

This example demonstrates how to isolate the emission of the Fe XVIII 93.9 Å line
in images from the AIA 94 Å channel.

The 94 Å channel contains a significant contribution from "warm" emission formed
at temperatures well below the roughly 7 MK at which Fe XVIII is formed.

Following the appendix of :cite:t:`warren2012`, we estimate the warm contribution
from a mixture of the 171 Å and 193 Å channels and subtract it from the 94 Å channel.
"""

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

import sunpy.map
from sunpy.data import manager

from aiapy.analysis import isolate_fe18_intensity
from aiapy.calibrate import correct_degradation, register, update_pointing
from aiapy.calibrate.utils import get_correction_table, get_pointing_table

###############################################################################
# The 94, 171, and 193 Å level 1 images are hosted in
# `aiapy-data <https://github.com/LM-SAL/aiapy-data>`__.
# SunPy's data manager downloads them once, verifies their SHA-256 hashes,
# and reuses the local cache on later runs.

base_url = "https://github.com/LM-SAL/aiapy-data/releases/download/v1/"


@manager.require(
    "aia_094_fe18",
    base_url + "aia.lev1_euv_12s.2011-11-08T163404Z.94.image_lev1.fits",
    "913aaa3b16fe74f524efffa40ab03df8a6d299243d974be32531c3b63f40e01e",
)
@manager.require(
    "aia_171_fe18",
    base_url + "aia.lev1_euv_12s.2011-11-08T163413Z.171.image_lev1.fits",
    "41810a805a30507d16417a8ccd36758a1a268671124912b7260264e91a6960d3",
)
@manager.require(
    "aia_193_fe18",
    base_url + "aia.lev1_euv_12s.2011-11-08T163409Z.193.image_lev1.fits",
    "715b42cb8097de0c0003e5ff79eaef16d0fcbb3817f2920844b6a0abf36088cd",
)
def load_images():
    """
    Load the three level 1 sample images.
    """
    return sunpy.map.Map([manager.get(name) for name in ("aia_094_fe18", "aia_171_fe18", "aia_193_fe18")])


map_94, map_171, map_193 = load_images()

###############################################################################
# We update pointing, register the images, correct degradation, and normalize
# by exposure time.

pointing_table = get_pointing_table("lmsal")
correction_table = get_correction_table(source="SSW")

level_1_5_maps = []
for level_1_map in (map_94, map_171, map_193):
    map_with_pointing = update_pointing(level_1_map, pointing_table=pointing_table)
    map_registered = register(map_with_pointing)
    map_corrected = correct_degradation(map_registered, correction_table=correction_table)
    level_1_5_maps.append(map_corrected / map_corrected.exposure_time)

###############################################################################
# The registered images can have slightly different dimensions, so we cut out
# the same region of the Sun from each of them.
#
# Figure 12 of :cite:t:`warren2012` shows this active region with a bright
# core of high-temperature emission.

bottom_left = SkyCoord(-285 * u.arcsec, 19 * u.arcsec, frame=level_1_5_maps[0].coordinate_frame)
top_right = SkyCoord(265 * u.arcsec, 569 * u.arcsec, frame=level_1_5_maps[0].coordinate_frame)
m_94, m_171, m_193 = [smap.submap(bottom_left, top_right=top_right) for smap in level_1_5_maps]

###############################################################################
# Now we can estimate and remove the warm contribution to the 94 Å channel.
# The difference between the observed 94 Å image and the estimated warm emission
# gives an estimate of the Fe XVIII 93.9 Å emission.

fe18_map = isolate_fe18_intensity(m_94, m_171, m_193)
warm_map = sunpy.map.Map(m_94.data - fe18_map.data, m_94.meta)

###############################################################################
# Let's compare the observed 94 Å image with the estimated warm emission and
# the isolated Fe XVIII emission.
#
# Subtracting the estimated warm emission highlights the compact core of
# high-temperature emission.
# All panels use the same intensity scale. Values outside the displayed range
# are clipped only for display; the maps retain their original values.

fig = plt.figure(figsize=(12, 4), layout="compressed")
fig.get_layout_engine().set(w_pad=0, wspace=0)
norm = ImageNormalize(m_94.data, vmin=0, interval=PercentileInterval(99.5), stretch=AsinhStretch(0.05), clip=True)
axes = []
for i, (smap, title) in enumerate(
    zip((m_94, warm_map, fe18_map), ("AIA 94 Å", "Warm 94 Å", "Fe XVIII 93.9 Å"), strict=True),
):
    ax = fig.add_subplot(1, 3, i + 1, projection=smap)
    image = smap.plot(axes=ax, annotate=False, norm=norm, cmap="sdoaia94")
    ax.set_title(title)
    ax.set_axis_off()
    axes.append(ax)

fig.colorbar(image, ax=axes, location="bottom", shrink=0.7, pad=0.03, label="Intensity [DN/s]", extend="both")
plt.show()

###############################################################################
# Note that this procedure cannot be applied to flare observations, where
# Fe XXIV 192.04 Å contributes to the 193 Å channel.
# The method is also unreliable in regions of very bright 1 MK emission,
# such as moss. The returned difference can be negative where the estimated
# warm contribution exceeds the observed 94 Å intensity.
# Averaging co-aligned images can improve the signal to noise; choose the
# averaging interval according to how quickly the source evolves.
