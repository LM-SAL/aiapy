"""
=============================================================
Deconvolving Images with the Instrument Point Spread Function
=============================================================

This example demonstrates how to deconvolve an AIA image with
the instrument point spread function (PSF).
"""
import matplotlib.pyplot as plt

import astropy.units as u
import sunpy.map
from astropy.coordinates import SkyCoord
from astropy.visualization import AsinhStretch, ImageNormalize, LogStretch
from sunpy.net import Fido, attrs

import aiapy.psf

#########################################
# AIA images are subject to convolution with the instrument point-spread
# function (PSF) due to effects introduced by the filter mesh of the telescope
# and the CCD, among others. This has the effect of "blurring" the image.
# The PSF diffraction pattern may also be particularly noticable during the
# impulsive phase of a flare where the intensity enhancement is very localized.
# To remove these artifacts, the PSF must be deconvolved from the image.
#
# First, we'll download a single level 1 image from the 171 Å channel from
# 7 June 2011. Note that deconvolution should be performed on level 1 images
# only. This is because, as with the level 1 data, the PSF model is defined
# on the CCD grid. Once deconvolved, the image can be passed to
# `~aiapy.calibrate.register`
# (see the :ref:`sphx_glr_generated_gallery_prepping_level_1_data.py` example).
q = Fido.search(
    attrs.Time('2011-06-07T06:52:00', '2011-06-07T06:52:10'),
    attrs.Instrument('AIA'),
    attrs.Wavelength(wavemin=171*u.angstrom, wavemax=171*u.angstrom),
)
m = sunpy.map.Map(Fido.fetch(q))
fig = plt.figure()
ax = fig.add_subplot(111, projection=m)
m.plot(axes=ax,)

#######################################
# Next, we'll calculate the PSF using `~aiapy.psf.psf` for the 171 Å channel.
# The PSF model accounts for several different effects, including diffraction
# from the mesh grating of the filters, charge spreading, and jitter. See
# `Grigis et al (2012) <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/psf/DOC/psfreport.pdf>`_
# for more details. Currently, this only works for
# :math:`4096\times4096` full frame images.
#
# Note that this will be significantly faster if you have a GPU and the `cupy`
# package installed.
psf = aiapy.psf.psf(m.wavelength)

#############################################
# We'll plot just a 500-by-500 pixel section centered on the center pixel. The
# diffraction "arms" extending from the center pixel can often be seen in
# flare observations due to the intense, small-scale brightening.
fov = 500
lc_x, lc_y = psf.shape[0]//2 - fov//2, psf.shape[1]//2 - fov//2
plt.imshow(psf[lc_x:lc_x+fov, lc_y:lc_y+fov],
           norm=ImageNormalize(vmin=1e-8, vmax=1e-3, stretch=LogStretch()))
plt.colorbar()
plt.show()

###############################################
# Now that we've downloaded our image and computed the PSF, we can deconvolve
# the image with the PSF using the
# `Richardson-Lucy deconvolution algorithm <https://en.wikipedia.org/wiki/Richardson%E2%80%93Lucy_deconvolution>`_.
# Note that passing in the PSF is optional. If you exclude it, it will be
# calculated automatically. However, when deconvolving many images of the same
# wavelength, it is most efficient to only calculate the PSF once.
#
# As with `~aiapy.psf.psf`, this will be much faster if you have
# a GPU and `cupy` installed.
m_deconvolved = aiapy.psf.deconvolve(m, psf=psf)

################################################
# Let's compare the convolved and deconvolved images.
norm = ImageNormalize(vmin=0, vmax=1.5e4, stretch=AsinhStretch(0.01))
fig = plt.figure()
ax = fig.add_subplot(121, projection=m)
m.plot(axes=ax, norm=norm)
ax = fig.add_subplot(122, projection=m_deconvolved)
m_deconvolved.plot(axes=ax, annotate=False, norm=norm)
ax.coords[0].set_axislabel(' ')
ax.coords[1].set_axislabel(' ')
ax.coords[1].set_ticklabel_visible(False)
plt.show()

#################################################
# The differences become a bit more obvious when we zoom in. Note that the
# deconvolution has the effect of "deblurring" the image.
left_corner = 500*u.arcsec, -600*u.arcsec
right_corner = 1000*u.arcsec, -100*u.arcsec
fig = plt.figure()
m_sub = m.submap(
    SkyCoord(*left_corner, frame=m.coordinate_frame),
    SkyCoord(*right_corner, frame=m.coordinate_frame)
)
ax = fig.add_subplot(121, projection=m_sub)
m_sub.plot(axes=ax, norm=norm)
m_deconvolved_sub = m_deconvolved.submap(
    SkyCoord(*left_corner, frame=m_deconvolved.coordinate_frame),
    SkyCoord(*right_corner, frame=m_deconvolved.coordinate_frame)
)
ax = fig.add_subplot(122, projection=m_deconvolved_sub)
m_deconvolved_sub.plot(axes=ax, annotate=False, norm=norm)
ax.coords[0].set_axislabel(' ')
ax.coords[1].set_axislabel(' ')
ax.coords[1].set_ticklabel_visible(False)
plt.show()
