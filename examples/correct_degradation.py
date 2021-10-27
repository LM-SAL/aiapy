"""
=====================================
Correcting for Instrument Degradation
=====================================

This example demonstrates the degradation of the filters on AIA over time.
"""

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u
from astropy.visualization import quantity_support, time_support
from sunpy.net import Fido
from sunpy.net import attrs as a

from aiapy.calibrate import degradation

# These are needed to allow the use of quantities and astropy
# time objects in the plot.
time_support(format='jyear')
quantity_support()

###########################################################
# The performance of the AIA telescope is unfortunately degrading over time,
# leading to the resulting images becoming increasingly dim. We
# can correct for this by modeling the degradation over time and
# then dividing the image intensity by this correction.
#
# First, let's fetch some metadata for the 335 Å channel of AIA between 2010
# and 2018 at a cadence of 30 days. We choose the 335 Å channel because it has experienced
# significant degradation compared to the other EUV channels.
results = Fido.search(
    a.Time('2010-06-01T00:00:00', '2021-06-01T00:00:00'),
    a.Sample(30*u.day),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.jsoc.Wavelength(335*u.angstrom),
)

###########################################################
# We only need the date and mean intensity columns from the
# metadata that was returned. We select those and nothing else.
table = results["jsoc"].show('DATE__OBS', 'DATAMEAN')
table['DATAMEAN'].unit = u.ct
table['DATE_OBS'] = astropy.time.Time(table['DATE__OBS'], scale='utc')
del table['DATE__OBS']

print(table)

###########################################################
# Next, we pass the date column to the `~aiapy.calibrate.correct_degradation`
# function. This function calculates the time-dependent correction factor
# based on the time and wavelength of the observation.
# We then divide the mean intensity by the correction factor to get the corrected intensity.
# For more details on how the correction factor is calculated, see the documentation for the
# `~aiapy.calibrate.degradation` function.
correction_factor = degradation(335*u.angstrom, table['DATE_OBS'])
# This correction can be applied to a sunpy map as well.
table['DATAMEAN_DEG'] = table['DATAMEAN'] / correction_factor

###########################################################
# To understand the effect of the degradation and the correction factor, we
# plot the corrected and uncorrected mean intensity as a function of time.
# Note that the uncorrected intensity decreases monotonically over time
# while the corrected intensity recovers to pre-2011 values in 2020.
plt.plot(table['DATE_OBS'], table['DATAMEAN'], label='mean', marker='o')
plt.plot(table['DATE_OBS'], table['DATAMEAN_DEG'], label='mean, corrected', marker='o')
plt.title(f'{(335*u.angstrom).to_string(format="latex")} Channel Degradation')
plt.legend(frameon=False)

plt.show()
