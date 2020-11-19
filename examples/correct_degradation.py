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
from astropy.visualization import time_support, quantity_support
from sunpy.net import attrs, jsoc

from aiapy.calibrate import degradation

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
channel = 335*u.angstrom
q = jsoc.JSOCClient().search(
    attrs.Time('2010-06-01T00:00:00', '2020-06-01T00:00:00'),
    attrs.Sample(30*u.day),
    jsoc.attrs.Series.aia_lev1_euv_12s,
    jsoc.attrs.Wavelength(channel),
)

###########################################################
# Select the date, mean, and max columns and clean up the table.
table = q.show('DATE__OBS', 'DATAMEAN')
table['DATAMEAN'].unit = u.ct
table['DATE_OBS'] = astropy.time.Time(table['DATE__OBS'], scale='utc')
del table['DATE__OBS']
print(table)

###########################################################
# Next, we can pass each map to the `~aiapy.calibrate.correct_degradation`
# function. This function calculates the time-dependent correction factor
# based on the time and wavelength of the observation, divides the intensity
# by the correction factor, and returns a new corrected map. For more details
# on how the correction factor is calculated, see the documentation for the
# `~aiapy.calibrate.degradation` function.
deg = degradation(channel, table['DATE_OBS'])
table['DATAMEAN_DEG'] = table['DATAMEAN'] / deg

###########################################################
# The effect of the degradation correction is more easily seen in the
# total flux as a function of time.
time_support(format='jyear')
quantity_support()
plt.plot(table['DATE_OBS'], table['DATAMEAN'], label='mean', marker='o')
plt.plot(table['DATE_OBS'], table['DATAMEAN_DEG'], label='mean, corrected', marker='o')
plt.title(f'{channel.to_string(format="latex")} Channel Degradation')
plt.legend(frameon=False)
plt.show()
