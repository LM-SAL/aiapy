.. _aiapy-prepping-level-1:

============================================
Preparing AIA data from level 1 to level 1.5
============================================

AIA data products provided by the JSOC are level 1 data products.
This means that the images still include the roll angle of the satellite and each channel may have a slightly different pixel scale.
Typically, before performing any sort of data analysis on AIA images, you will want to promote your AIA data from level 1 to level 1.5.
This is important if you want to compare images from different channels or create Differential Emission Measure (DEM) maps.

The promotion to level 1.5 involves updating the pointing keywords, removing the roll angle, scaling the image to a resolution of 0.6 arcseconds per pixel, and translating the image such that the center of the Sun is located in the center of the image.

In IDL, this is done with the ``aia_prep.pro`` procedure in SSWIDL as described in the `SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`__.
The following example, :ref:`sphx_glr_generated_gallery_prepping_level_1_data.py` demonstrates how to achieve this in Python with ``aiapy``.

There are also additional processing steps that can be applied to the level 1 AIA images.
If you want to do any additional data processing steps (e.g., PSF deconvolution) should be done in the following order:

1. Pointing correction (`aiapy.calibrate.update_pointing`)
2. Image respiking (`aiapy.calibrate.respike`)
3. PSF deconvolution (`aiapy.psf.deconvolve`)
4. Registration (`aiapy.calibrate.register`)
5. Degradation correction (`aiapy.calibrate.correct_degradation`)
6. Exposure normalization

.. note::

   * Level 1.5, in its typical usage, only includes steps 1 and 4.
     Unless stated otherwise, any science publication mentioning level 1.5 AIA data does not include steps 2, 3, 5 and 6.
   * The PSF functions are defined on the level 1 pixel grid so PSF deconvolution **MUST** be done on the level 1 data products (i.e., before image registration).
     This is described in the PSF gallery example :ref:`sphx_glr_generated_gallery_skip_psf_deconvolution.py`.
   * The pointing update should be done prior to image registration as the updated keywords, namely ``CRPIX1`` and ``CRPIX2``, are used in the image registration step.
   * The exposure time normalization and degradation correction (`aiapy.calibrate.correct_degradation`) operations are just scalar multiplication and are thus linear such that their ordering is inconsequential.
   * Exposure time normalization can be performed by simply dividing a map by the exposure time property, ``my_map / my_map.exposure_time``.
