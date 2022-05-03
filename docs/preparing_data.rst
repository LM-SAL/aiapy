Preparing AIA data
==================

The common usecase for ``aiapy`` is to transform level 1 AIA data to level 1.5.
This is called `aia_prep.pro` procedure in SSWIDL as described in the `SDO Analysis Guide <https://www.lmsal.com/sdodocs/doc/dcur/SDOD0060.zip/zip/entry/index.html>`__.

The following example, :ref:`sphx_glr_generated_gallery_prepping_level_1_data.py` showcases how to mimic this in Python via ``aiapy``.

If you want to do more advanced things, for example, a PSF deconvolution.
You will want to do it in the following order:

1. PSF deconvolution (`aiapy.psf.deconvolve`)
2. Pointing correction (`aiapy.calibrate.update_pointing`)
3. Registration (`aiapy.calibrate.register`)
4. Degradation correction (`aiapy.calibrate.correct_degradation`)
5. Exposure normalization (`aiapy.calibrate.normalize_exposure`)

A few notes on this:

* The PSF functions are defined on the level 1 pixel grid so PSF deconvolution MUST be done on the level 1 data products (i.e. before image registration).
  This is described in the PSF gallery example :ref:`sphx_glr_generated_gallery_skip_psf_deconvolution.py`.
* The pointing update should be done prior to image registration as the updated keywords, namely ``CRPIX1`` and ``CRPIX2``, are used in the image registration step.
  More details can be found in this gallery example :ref:`sphx_glr_generated_gallery_update_header_keywords.py`.
* The exposure time normalization (`aiapy.calibrate.normalize_exposure`) and degradation correction (`aiapy.calibrate.correct_degradation`) operations are just scalar multiplication and are thus linear such that their ordering is inconsequential.
* Level 1.5, in its typical usage, only includes steps 2 and 3.
  Unless stated otherwise, a science publication mentioning level 1.5 AIA data does not include steps 1, 4 and 5.
