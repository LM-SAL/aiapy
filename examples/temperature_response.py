"""
========================================
Computing temperature response functions
========================================

This example shows how to compute the temperature response functions of the
AIA EUV channels, as SolarSoft's ``aia_get_response(/temperature, /dn)`` does.
"""

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u

from aiapy.response import Channel, get_source_spectra, temperature_response

###############################################################################
# The temperature response folds the emissivity of an optically thin plasma
# with the wavelength response of a channel. `aiapy.response.get_source_spectra`
# returns the emissivity that SolarSoft uses, which is downloaded the first time.

spectra = get_source_spectra()
print(spectra.meta["chianti_version"], spectra.meta["abundance"])

###############################################################################
# These are the responses of the six channels usually used for DEM analysis,
# normalized to EVE and corrected for the degradation at a given time.

obstime = astropy.time.Time("2024-05-10T00:00:00", scale="utc")

fig, ax = plt.subplots()
for wavelength in [94, 131, 171, 193, 211, 335] * u.angstrom:
    channel = Channel(wavelength, include_eve_correction=True)
    temperature, response = temperature_response(channel, spectra, obstime)
    ax.plot(temperature, response, label=f"{wavelength.value:.0f} Å")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e5, 1e8)
ax.set_ylim(1e-28, 1e-23)
ax.set_xlabel(f"$T$ [{temperature.unit.to_string('latex')}]")
ax.set_ylabel(f"$R(T)$ [{response.unit.to_string('latex')}]")
ax.legend(frameon=False)

###############################################################################
# SolarSoft also has an empirical correction to the 94 Å response for emission
# missing from CHIANTI, the "chiantifix".

channel = Channel(94 * u.angstrom, include_eve_correction=True)
temperature, response = temperature_response(channel, spectra, obstime)
_, response_chiantifix = temperature_response(channel, spectra, obstime, chiantifix=True)

fig, ax = plt.subplots()
ax.plot(temperature, response, label="94 Å")
ax.plot(temperature, response_chiantifix, label="94 Å with chiantifix")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e5, 1e8)
ax.set_ylim(1e-30, 1e-26)
ax.set_xlabel(f"$T$ [{temperature.unit.to_string('latex')}]")
ax.set_ylabel(f"$R(T)$ [{response.unit.to_string('latex')}]")
ax.legend(frameon=False)

###############################################################################
# An array of times gives one response per time, here showing how the
# 335 Å channel has degraded.

obstimes = astropy.time.Time([f"{year}-01-01T00:00:00" for year in range(2011, 2026, 2)], scale="utc")
channel = Channel(335 * u.angstrom, include_eve_correction=True)
temperature, response = temperature_response(channel, spectra, obstimes)

fig, ax = plt.subplots()
for t, row in zip(obstimes, response, strict=True):
    ax.plot(temperature, row, label=t.strftime("%Y"))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e5, 1e8)
ax.set_ylim(1e-29, 1e-26)
ax.set_xlabel(f"$T$ [{temperature.unit.to_string('latex')}]")
ax.set_ylabel(f"$R(T)$ [{response.unit.to_string('latex')}]")
ax.legend(frameon=False)

plt.show()
