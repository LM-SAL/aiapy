"""
===================================================
Isolating Fe XVIII emission in the AIA 94 Å channel
===================================================

The AIA 94 Å channel mixes Fe XVIII 93.9 Å emission, formed near 7 MK, with
cooler emission. Following :cite:t:`warren2012`, this example estimates and
subtracts the warm contribution using 171 Å and 193 Å images.

.. warning::

    These fixed coefficients were fitted to 2010 March 22 observations. They
    may not suit other dates or calibrations; check the warm estimate for your data.
"""

import matplotlib.pyplot as plt
import pooch

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

import sunpy.map

from aiapy.analysis import isolate_fe18_intensity
from aiapy.calibrate import correct_degradation, register, update_pointing
from aiapy.calibrate.utils import get_correction_table, get_pointing_table

###############################################################################
# We require some SDO/AIA data, for this examples we have hosted some sample data
# but you will want to use SunPy's Fido.

base_url = "https://github.com/LM-SAL/aiapy-data/releases/download/v1/"
file_94 = pooch.retrieve(
    base_url + "aia.lev1_euv_12s.2011-11-08T163404Z.94.image_lev1.fits",
    known_hash="913aaa3b16fe74f524efffa40ab03df8a6d299243d974be32531c3b63f40e01e",
)
file_171 = pooch.retrieve(
    base_url + "aia.lev1_euv_12s.2011-11-08T163413Z.171.image_lev1.fits",
    known_hash="41810a805a30507d16417a8ccd36758a1a268671124912b7260264e91a6960d3",
)
file_193 = pooch.retrieve(
    base_url + "aia.lev1_euv_12s.2011-11-08T163409Z.193.image_lev1.fits",
    known_hash="715b42cb8097de0c0003e5ff79eaef16d0fcbb3817f2920844b6a0abf36088cd",
)

map_94, map_171, map_193 = sunpy.map.Map([file_94, file_171, file_193])

###############################################################################
# Update pointing, register, correct degradation, and divide by exposure time.

pointing_table = get_pointing_table("lmsal")
correction_table = get_correction_table(source="SSW")

level_1_5_maps = []
for level_1_map in (map_94, map_171, map_193):
    map_with_pointing = update_pointing(level_1_map, pointing_table=pointing_table)
    map_registered = register(map_with_pointing)
    map_corrected = correct_degradation(map_registered, correction_table=correction_table)
    level_1_5_maps.append(map_corrected / map_corrected.exposure_time)

###############################################################################
# Registered images can have different shapes, so crop to the same solar region.
#
# The function does not check for alignment or different observation timing.

bottom_left = SkyCoord(-285 * u.arcsec, 19 * u.arcsec, frame=level_1_5_maps[0].coordinate_frame)
top_right = SkyCoord(265 * u.arcsec, 569 * u.arcsec, frame=level_1_5_maps[0].coordinate_frame)
m_94, m_171, m_193 = [smap.submap(bottom_left, top_right=top_right) for smap in level_1_5_maps]

###############################################################################
# Subtract the warm estimate to isolate Fe XVIII emission.

fe18_map = isolate_fe18_intensity(m_94, m_171, m_193)
warm_map = sunpy.map.Map(m_94.data - fe18_map.data, m_94.meta, mask=fe18_map.mask)

###############################################################################
# Compare the observed 94 Å, estimated warm, and Fe XVIII images.
#
# The subtraction highlights the compact hot core.

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
# During flares, Fe XXIV 192.04 Å can enter the 193 Å channel and bright moss
# makes the estimate unreliable. The residual can be negative where the warm
# estimate exceeds the observed 94 Å intensity.
# One can average aligned images to reduce noise.
