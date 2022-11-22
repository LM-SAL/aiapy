r"""
=================================
Converting DN/s to Physical Units
=================================

This example demonstrates how to convert AIA Level 1 images that are in DN/s to physical units, A.

One might ask why not a function, the reason being is that several assumptions have to be made in.
As a result, this example will detail each assumption, the reason for it.
This will allow a more informed user to change each assumption made to suit their own workflows.

This example will make heavy use of the equations found in "Initial Calibration of the Atmospheric Imaging Assembly (AIA) on the Solar Dynamics Observatory (SDO)" by `Boerner et al. (2012) <https://link.springer.com/article/10.1007/s11207-011-9804-8>`__.


:math:`P_{c}=\int d\lambda R_{c}(\lambda) I(\lambda)`

:math:`P_{c} \approx (FWHM)_{\lambda} \tilde{R}_{c} \tilde{I} \therefore \tilde{I} \approx P_{c} / (FWHM)_{\lambda} \tilde{R}_{c}`

:math:`\frac{[DN \, pix^{-1} \, s^{-1}]}{[\AA][DN \, ph^{-1} \, sr \, pix^{-1}]} \Rightarrow \frac{ph}{s \, \AA \, sr}`

"""
import matplotlib.pyplot as plt

import astropy.units as u
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import aiapy.data.sample as sample_data
from aiapy.calibrate import normalize_exposure, register, update_pointing

###########################################################
# To perform unit conversation, we will start with a level 1 FITS file.
# We will use a level 1 171 Å AIA image from the `aiapy` sample data and
# read it using a `~sunpy.map.Map`.

aia_map = sunpy.map.Map(sample_data.AIA_171_IMAGE)

###########################################################
# We will be doing operations on the data array of the map.
# However, we want to keep track of the metadata and the units.

print(aia_map.unit)
# We can access a unit-ed array calling: aia_map.quantity

###########################################################
# The unit on the map is ``u.ct`` which is counts and not digital number [DN].
# This might be confusing but exists for two reasons.
#
# 1) `astropy` did not add DN as a unit until version 4.3 (2021).
# 2) It is not a valid FITS unit and as such, you can not save
#    out a FITS file with DN as a unit.
#
# Due to this, counts was chosen as the unit by default.
# However, the data is in DN and not counts.
#
# For now we will normalize to the exposure time to get to "DN"/s.

normalized_map = normalize_exposure(aia_map)
print(normalized_map.unit)

###########################################################
# Let us now begin the process of converting DN to physical units.
#
# "Initial Calibration of the Atmospheric Imaging Assembly (AIA) on the Solar Dynamics Observatory (SDO)"
# by `Boerner et al. (2012) <https://link.springer.com/article/10.1007/s11207-011-9804-8>`__
# Will be the basis for this example.
#
# Starting with Equation 1:
#
# :math:`p_{i}(\bold{x}) = \int_{0}^{\infty} \eta_{i}(\lambda)d\lambda \int_{pixel\,\bold{x}} I(\lambda, \theta)d\theta`
#
# Here :math:`\eta_{i}` is the efficiency function of the :math:`i`th channel of the telescope in units
# of DN per unit flux at the aperture
#
# # :math: p_{i}(\bold{x}) = \int_{0}^{\infty} \frac{I_{i}(\bold{x}, \lambda) \phi(\lambda) d\lambda}{h \nu_{i}(\bold{x}) c}
# :math:`P_{c}=\int d\lambda R_{c}(\lambda) I(\lambda)`

# map.quantity / (channel.wavelenth_response * channel.plate_scale)
# Wavelength degedration?
# FWHM astroy model
# Get models to fit as well?!
# Check how xrtpy do temperature responses
# Add photosphere response files
