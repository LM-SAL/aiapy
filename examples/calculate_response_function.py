"""
=======================================
Computing wavelength response functions
=======================================

This example shows how to compute the
wavelength response function of the 335 Å channel as
well as explore the different properties of the
telescope channels.
"""

import matplotlib.pyplot as plt

import astropy.time
import astropy.units as u

from aiapy.calibrate.util import get_correction_table
from aiapy.response import Channel

###############################################################################
# Since AIA uses narrow-band filters, other wavelengths (outside of the nominal
# wavelength attributed to each filter) contribute to the image data.
# Computing these response functions allow us to see which other wavelengths
# contribute to the total intensity in each image.
#
# First, create a `aiapy.response.Channel` object by specifying the
# wavelength of the channel. In this case, we'll
# choose the 335 Å channel, but this same workflow
# can be applied to any of the EUV or UV channels
# on AIA. This may take a few seconds the first time you do
# this as the most recent instrument data file will
# need to be downloaded from a remote server. Subsequent
# calls will know that the data has been downloaded.

aia_335_channel = Channel(335 * u.angstrom)

###############################################################################
# From `Boerner et al. (2012) <https://doi.org/10.1007/s11207-011-9804-8>`_,
# the wavelength response function is given by,
#
# .. math::
#   R(\lambda) = A_{geo}R_P(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)
#   D(\lambda)Q(\lambda)G(\lambda),
#
# where
#
# - :math:`A_{geo}` geometrical collecting area
# - :math:`R_P,R_S` reflectances of primary and secondary mirrors, respectively
# - :math:`T_E, T_F` transmission efficiency of the entrance and focal-plane
#   filters, respectively
# - :math:`D` contaminant transmittance of optics
# - :math:`Q` quantum efficiency of the CCD
# - :math:`G` gain of the CCD camera system
#
# The `aiapy.response.Channel` object provides an interface to all of these
# properties of the telescope. Below, we show how to plot several of these
# properties as a function of wavelength.

# Reflectance
fig = plt.figure()
ax = fig.add_subplot(221)
ax.plot(aia_335_channel.wavelength, aia_335_channel.primary_reflectance, label=r"$R_P$")
ax.plot(aia_335_channel.wavelength, aia_335_channel.secondary_reflectance, label=r"$R_S$")
ax.set_ylabel(r"Reflectance")
ax.set_xlim(50, 400)
ax.set_xlabel(r"$\lambda$ [Å]")
ax.legend(frameon=False)

# Transmittance
ax = fig.add_subplot(222)
ax.plot(aia_335_channel.wavelength, aia_335_channel.entrance_filter_efficiency, label=r"$T_E$")
ax.plot(aia_335_channel.wavelength, aia_335_channel.focal_plane_filter_efficiency, label=r"$T_F$")
ax.set_ylabel(r"Transmittance")
ax.set_xlim(50, 400)
ax.set_xlabel(r"$\lambda$ [Å]")
ax.legend(frameon=False)

# Contamination
ax = fig.add_subplot(223)
ax.plot(aia_335_channel.wavelength, aia_335_channel.contamination)
ax.set_ylabel(r"Contamination, $D(\lambda)$")
ax.set_xlim(50, 400)
ax.set_xlabel(r"$\lambda$ [Å]")

# Quantumn efficiency
ax = fig.add_subplot(224)
ax.plot(aia_335_channel.wavelength, aia_335_channel.quantum_efficiency)
ax.set_ylabel(r"Quantum Efficiency, $Q(\lambda)$")
ax.set_xlim(50, 800)
ax.set_xlabel(r"$\lambda$ [Å]")

plt.tight_layout()

###############################################################################
# Additionally, `aiapy.response.Channel` provides a method for calculating
# the wavelength response function using the equation above,

correction_table = get_correction_table("jsoc")
wavelength_response_335 = aia_335_channel.wavelength_response(correction_table=correction_table)
print(wavelength_response_335)

###############################################################################
# We can then plot the response as a function of
# wavelength.

fig = plt.figure()

ax = fig.gca()
ax.plot(aia_335_channel.wavelength, wavelength_response_335)
ax.set_xlim((aia_335_channel.channel + [-10, 10] * u.angstrom).value)
ax.set_ylim(0, 0.03)
ax.set_xlabel(r"$\lambda$ [Å]")
ax.set_ylabel(f'$R(\\lambda)$ [{wavelength_response_335.unit.to_string("latex")}]')

###############################################################################
# On telescopes 1, 3, and 4, both channels are always illuminated.
# This can lead to "crosstalk" contamination in a channel from the channel with
# which it shares a telescope. This impacts the 94 Å and 304 Å channels
# as well as the 131 Å and 335 Å channels. This effect is included
# by default in the wavelength response calculation. To exclude this
# effect,

wavelength_response_335_no_cross = aia_335_channel.wavelength_response(
    include_crosstalk=False, correction_table=correction_table
)

###############################################################################
# If we look at the response around 131 Å (the channel with which 335 Å shares
# a telescope), we can see the effect that the channel crosstalk has on the
# 335 Å response function.

fig = plt.figure()

ax = fig.gca()
ax.plot(aia_335_channel.wavelength, wavelength_response_335, label="crosstalk")
ax.plot(aia_335_channel.wavelength, wavelength_response_335_no_cross, label="no crosstalk")
ax.set_xlim(50, 350)
ax.set_xlabel(r"$\lambda$ [Å]")
ax.set_ylabel(f'$R(\\lambda)$ [{wavelength_response_335.unit.to_string("latex")}]')
ax.legend(loc=1, frameon=False)

###############################################################################
# We can also incorporate various corrections to the
# response functions, including a time-dependent
# degradation correction as well as a correction based
# on the EVE calibration. The latter also includes the
# time-dependent correction. As an example, to apply the
# two aforementioned corrections given the degradation as
# of 1 January 2019,

obstime = astropy.time.Time("2019-01-01T00:00:00")
wavelength_response_335_time = aia_335_channel.wavelength_response(obstime=obstime, correction_table=correction_table)
wavelength_response_335_eve = aia_335_channel.wavelength_response(
    obstime=obstime, include_eve_correction=True, correction_table=correction_table
)

###############################################################################
# We can then compare the two corrected response
# functions to the uncorrected case.

fig = plt.figure()
ax = fig.gca()
ax.plot(aia_335_channel.wavelength, wavelength_response_335, label="uncorrected")
ax.plot(aia_335_channel.wavelength, wavelength_response_335_time, label="degradation correction")
ax.plot(aia_335_channel.wavelength, wavelength_response_335_eve, label="EVE correction")
ax.set_xlim((aia_335_channel.channel + [-20, 20] * u.angstrom).value)
ax.set_ylim(0, 0.03)
ax.set_xlabel(r"$\lambda$ [Å]")
ax.set_ylabel(f'$R(\\lambda)$ [{wavelength_response_335.unit.to_string("latex")}]')
ax.legend(loc=2, frameon=False)

plt.show()
