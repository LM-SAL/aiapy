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
# First, let's fetch one 335 Å AIA observational metadata between 2010
# and 2018, one every 90 days.
# We choose the 335 Å channel because it has experienced
# significant degradation compared to the other EUV channels.
results = Fido.search(
    a.Time('2010-06-01T00:00:00', '2021-06-01T00:00:00'),
    a.Sample(90*u.day),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.jsoc.Wavelength(335*u.angstrom),
)

###########################################################
# We only need some of the metadata that was returned.
# In this case, the date, mean, and max columns and nothing else.
table = results["jsoc"].show('DATE__OBS', 'DATAMEAN')
table['DATAMEAN'].unit = u.ct
table['DATE_OBS'] = astropy.time.Time(table['DATE__OBS'], scale='utc')
del table['DATE__OBS']

print(table)

###########################################################
# Next, we can pass the date column to the `~aiapy.calibrate.correct_degradation`
# function. Which calculates the time-dependent correction factor
# based on the time and wavelength of the observation, divides the intensity
# by the correction factor, and returns correction factor. For more details
# on how the correction factor is calculated, see the documentation for the
# `~aiapy.calibrate.degradation` function.
correction_factor = degradation(335*u.angstrom, table['DATE_OBS'])
# This correction can be applied to a sunpy map as well.
table['DATAMEAN_DEG'] = table['DATAMEAN'] / correction_factor

###########################################################
# To showcase the effect of degradation, using the total
# flux as a function of time is the best way to visualize it.
plt.plot(table['DATE_OBS'], table['DATAMEAN'], label='mean', marker='o')
plt.plot(table['DATE_OBS'], table['DATAMEAN_DEG'], label='mean, corrected', marker='o')
plt.title(f'{(335*u.angstrom).to_string(format="latex")} Channel Degradation')
plt.legend(frameon=False)

plt.show()
