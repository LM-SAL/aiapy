"""
=====================================
Correcting for instrument degradation
=====================================

This example demonstrates the degradation of the filters on AIA and how to correct it.
"""

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u
from astropy.visualization import time_support

from sunpy.net import Fido
from sunpy.net import attrs as a

from aiapy.calibrate import degradation
from aiapy.calibrate.util import get_correction_table

# This lets you pass `astropy.time.Time` objects directly to matplotlib
time_support(format="jyear")

###############################################################################
# First, let's fetch the metadata for the 335 Å channel of AIA between 2021
# and 2023 at a cadence of 7 days. We choose the 335 Å channel because it has experienced
# significant degradation compared to the other EUV channels.

results = Fido.search(
    a.Time("2021-01-01T00:00:00", "2023-01-01T00:00:00"),
    a.Sample(7 * u.day),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.jsoc.Wavelength(335 * u.angstrom),
)

###############################################################################
# We only need the date and mean intensity columns from the
# metadata that was returned. We select those and nothing else.

table = results["jsoc"].show("DATE__OBS", "DATAMEAN")
table["DATAMEAN"].unit = u.DN
table["DATE_OBS"] = astropy.time.Time(table["DATE__OBS"], scale="utc")
del table["DATE__OBS"]

print(table)

###############################################################################
# Next, we pass the date column to the `aiapy.calibrate.correct_degradation`
# function. This function calculates the time-dependent correction factor
# based on the time and wavelength of the observation.
# We then divide the mean intensity by the correction factor to get the corrected intensity.
# For more details on how the correction factor is calculated, see the documentation for the
# `aiapy.calibrate.degradation` function.

correction_factor = degradation(335 * u.angstrom, table["DATE_OBS"], correction_table=get_correction_table("jsoc"))
table["DATAMEAN_DEG"] = table["DATAMEAN"] / correction_factor

###############################################################################
# To understand the effect of the degradation and the correction factor, we
# plot the corrected and uncorrected mean intensity as a function of time.
# Note that the uncorrected intensity decreases monotonically over time
# while the corrected intensity recovers to pre-2011 values in 2020.

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(table["DATE_OBS"], table["DATAMEAN"], label="mean", marker="o")
ax.plot(table["DATE_OBS"], table["DATAMEAN_DEG"], label="corrected mean", marker="+")
ax.set_title(f'{(335*u.angstrom).to_string(format="latex")} Channel Degradation')
ax.legend(frameon=False)

plt.show()
