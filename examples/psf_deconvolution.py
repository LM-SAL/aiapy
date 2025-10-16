"""
===================================================================
Deconvolving images with the instrument Point Spread Function (PSF)
===================================================================

This example demonstrates how to deconvolve an AIA image with
the instrument point spread function (PSF).
"""
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, LogStretch

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time

import aiapy.psf

###############################################################################
# AIA images are subject to convolution with the instrument point-spread
# function (PSF) due to effects introduced by the filter mesh of the telescope
# and the CCD, among others. This has the effect of "blurring" the image.
# The PSF diffraction pattern may also be particularly noticeable during the
# impulsive phase of a flare where the intensity enhancement is very localized.
# To remove these artifacts, the PSF must be de-convolved from the image.
#
# First, we'll use a single level 1 image from the 171 Å channel from
# 10 September 2017. Note that deconvolution should be performed on level 1 images
# only. This is because, as with the level 1 data, the PSF model is defined
# on the CCD grid. Once deconvolved, the image can be passed to
# `aiapy.calibrate.register`
# (see the :ref:`sphx_glr_generated_gallery_prepping_level_1_data.py` example).

t_start = parse_time("2017-09-10T20:00:00")
search_results = Fido.search(
    a.Time(t_start, t_start + 12 * u.s),
    a.Instrument.aia,
    a.Wavelength(171 * u.angstrom),
)
files = Fido.fetch(search_results)
aia_map = sunpy.map.Map(files[0])

###############################################################################
# Next, we'll calculate the PSF using `aiapy.psf.calculate_psf` for the 171 Å channel.
# The PSF model accounts for several different effects, including diffraction
# from the mesh grating of the filters, charge spreading, and jitter. See
# `Grigis et al (2012) <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/psf/DOC/psfreport.pdf>`_
# for more details. Currently, this only works for  :math:`4096\times4096` full frame images.
#
# Note that this will be significantly faster if you have installed `jax <https://docs.jax.dev/en/latest/>`__.

psf = aiapy.psf.calculate_psf(aia_map.wavelength)

###############################################################################
# We'll plot just a 500-by-500 pixel section centered on the center pixel. The
# diffraction "arms" extending from the center pixel can often be seen in
# flare observations due to the intense, small-scale brightening.

fov = 500
lc_x, lc_y = psf.shape[0] // 2 - fov // 2, psf.shape[1] // 2 - fov // 2
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.imshow(
    psf[lc_x : lc_x + fov, lc_y : lc_y + fov],
    norm=ImageNormalize(vmin=1e-8, vmax=1e-3, stretch=LogStretch()),
    origin="lower",
)
ax.set_title(f"AIA {aia_map.wavelength.value:.0f} Å PSF")
ax.set_xlabel("Pixels")
ax.set_ylabel("Pixels")
plt.colorbar(ax.images[0], ax=ax, label="Normalized Intensity")

###############################################################################
# Now that we've downloaded our image and computed the PSF, we can deconvolve
# the image with the PSF using the
# `Richardson-Lucy deconvolution algorithm <https://en.wikipedia.org/wiki/Richardson%E2%80%93Lucy_deconvolution>`_.
# Note that passing in the PSF is optional. If you exclude it, it will be
# calculated automatically. However, when deconvolving many images of the same
# wavelength, it is most efficient to only calculate the PSF once.
#
# As with `aiapy.psf.calculate_psf`, this will be much faster if you have
# a Nvidia GPU and `cupy` installed.

aia_map_deconvolved = aiapy.psf.deconvolve(aia_map, psf=psf)

###############################################################################
# Let's compare the convolved and deconvolved images. We will now zoom in on some
# loop structures that extend off the limb to examine the differences between
# our image before and after deconvolution.
#
# The differences become a bit more obvious when we zoom in. Note that the
# deconvolution has the effect of "deblurring" the image.

bottom_left = SkyCoord(750, -375, unit="arcsec", frame=aia_map.coordinate_frame)
fov = {"width": 400 * u.arcsec, "height": 400 * u.arcsec}

aia_map_sub = aia_map.submap(bottom_left, **fov)
aia_map_deconvolved_sub = aia_map_deconvolved.submap(bottom_left, **fov)

fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(121, projection=aia_map_sub)
aia_map_sub.plot(axes=ax1, title="Before Deconvolution")
ax2 = fig.add_subplot(122, projection=aia_map_deconvolved_sub)
aia_map_deconvolved_sub.plot(axes=ax2, title="After Deconvolution")
# Just remove the ticks and labels
lon = ax2.coords[0]
lat = ax2.coords[1]
lon.set_ticks_visible(False)
lon.set_ticklabel_visible(False)
lat.set_ticks_visible(False)
lat.set_ticklabel_visible(False)
lon.set_axislabel("")
lat.set_axislabel("")
fig.tight_layout()

###############################################################################
# Note that comparing the images on the left and right, deconvolving
# with the PSF removes the appearance of the diffraction pattern near the
# bright loop top. This also brings out more detail in the strands of the
# post-flare arcade.
#
# We can see this more clearly by plotting the intensity of the original
# and deconvolved images across the top of the loop top.

line_coords = SkyCoord([950, 1050], [-137, -137], unit=(u.arcsec, u.arcsec), frame=aia_map_sub.coordinate_frame)
intensity_coords = sunpy.map.pixelate_coord_path(aia_map_sub, line_coords)
before_intensity = sunpy.map.sample_at_coords(aia_map_sub, intensity_coords)
after_intensity = sunpy.map.sample_at_coords(aia_map_deconvolved_sub, intensity_coords)
angular_separation = intensity_coords.separation(intensity_coords[0]).to(u.arcsec)

fig = plt.figure(figsize=(12, 5))
ax1 = fig.add_subplot(121, projection=aia_map_sub)
aia_map_sub.plot(axes=ax1)
ax1.plot_coord(intensity_coords)

ax2 = fig.add_subplot(122)
ax2.plot(angular_separation, before_intensity, label="Original")
ax2.plot(angular_separation, after_intensity, label="Deconvolved")
ax2.set_xlabel("Angular distance along slit [arcsec]")
ax2.set_ylabel(f"Intensity [{aia_map_sub.unit}]")
ax2.legend(loc="best", ncol=2, frameon=False)
fig.tight_layout()

plt.show()
