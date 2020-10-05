"""
=====================================
Correcting for Instrument Degradation
=====================================

This example demonstrates how to correct an AIA image
to account for the degradation of the telescopes over time.
"""

import matplotlib.pyplot as plt
import astropy.units as u
import astropy.time
from astropy.visualization import ImageNormalize, SqrtStretch, time_support
from sunpy.net import Fido, attrs
import sunpy.map

from aiapy.calibrate import correct_degradation
from aiapy.calibrate.util import get_correction_table


###########################################################
# The performance of the AIA telescope is degrading over time,
# leading to the resulting images becoming increasingly dim. We
# can correct for this by modeling the degradation over time and
# then dividing the image intensity by this correction.
#
# First, let's fetch one 335 Å AIA observation from the
# `Virtual Solar Observatory <https://sdac.virtualsolar.org/cgi/search>`_
# for every other year between 2010 and 2018 and create a list of
# `~sunpy.map.Map` objects. We choose the 335 Å channel because it has
# experienced significant degradation compared to the other EUV channels.
q = Fido.search(
    attrs.Time('2010-06-01T00:00:00', '2018-06-01T00:00:00'),
    attrs.Sample(2*u.year),
    attrs.Instrument('AIA'),
    attrs.Wavelength(335*u.angstrom),
)
maps = sunpy.map.Map(sorted(Fido.fetch(q)))

###########################################################
# Next, we can pass each map to the `~aiapy.calibrate.correct_degradation`
# function. This function calculates the time-dependent correction factor
# based on the time and wavelength of the observation, divides the intensity
# by the correction factor, and returns a new corrected map. For more details
# on how the correction factor is calculated, see the documentation for the
# `~aiapy.calibrate.degradation` function.
correction_table = get_correction_table()
maps_corrected = [correct_degradation(m, correction_table=correction_table) for m in maps]

###########################################################
# Let's plot the uncorrected and corrected images for each year to show the
# degradation over time.
norm = ImageNormalize(vmin=0, vmax=1e2, stretch=SqrtStretch())
fig = plt.figure(figsize=(len(maps), 2))
for i, (m, mc) in enumerate(zip(maps, maps_corrected)):
    ax = fig.add_subplot(2, len(maps), i+1, projection=m)
    m.plot(axes=ax, norm=norm, annotate=False)
    ax.set_title(m.date.datetime.year)
    ax.coords[0].set_ticks_visible(False)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticks_visible(False)
    ax.coords[1].set_ticklabel_visible(False)
    ax = fig.add_subplot(2, len(maps), i+1+len(maps), projection=mc)
    mc.plot(axes=ax, norm=norm, annotate=False,)
    ax.coords[0].set_ticks_visible(False)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticks_visible(False)
    ax.coords[1].set_ticklabel_visible(False)

###########################################################
# The effect of the degradation correction is more easily seen in the
# total flux as a function of time.
flux_corrected = u.Quantity([m.data.sum() for m in maps_corrected])
flux_uncorrected = u.Quantity([m.data.sum() for m in maps])
time = astropy.time.Time([m.date for m in maps])
time_support()
fig = plt.figure()
ax = fig.gca()
ax.plot(time, flux_uncorrected, label='uncorrected', marker='o')
ax.plot(time, flux_corrected, label='corrected', marker='o')
ax.set_xlabel('Time')
ax.set_ylabel('Total Intensity [DN]')
ax.legend(frameon=False)
