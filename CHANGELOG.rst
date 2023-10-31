Aiapy 0.7.4 (2023-10-31)
========================

Bug Fixes
---------

- Fixed mismatch with the sample data downloder.
- Fixed theme build to use the new sunpy theme.


0.7.3 (2023-04-05)
==================

Bug Fixes
---------

- Fixed missing citation. (`#177 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/177>`__)


0.7.2 (2023-01-18)
==================

Breaking Changes
----------------

- Removed kwargs from ``correct_degradation``, ``degradation`` to ensure that the correct keywords are passed into these functions and the function calls that require these keywords. (`#170 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/170>`__)
- Several private internal functions now raise `ValueError` instead of `IndexError`, some of these will now be raised to the user when calling ``correct_degradation``, therefore any code that checks for this error type will need updating.  (`#170 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/170>`__)

Bug Fixes
---------

- Improve error message for degradation with incorrect response table version. (`#170 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/170>`__)
- Fixed download issue for error tables and parfive due to a mirror change. (`#172 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/172>`__)

Documentation
-------------

- Added a warning about ``register`` that it will not return a 4096 by 4096 Map, but most likely a 4094 by 4094 Map. (`#170 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/170>`__)

Trivial Changes
---------------

- Added several tests to ensure that ``degradation`` works on all wavelengths and the ones it does not, raise the correct error. (`#170 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/170>`__)


0.7.1 (2022-11-28)
==================

Bug Fixes
---------

- Change SSW mirror due to the old one being down. (`#167 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/167>`__)


0.7 (2022-08-04)
================

Breaking Changes
----------------

- Minimum version of ``Python`` supported is now ``Python`` 3.8.0. (`#159 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/159>`__)
- Minimum version of ``sunpy`` supported is now ``sunpy`` 4.0.0 LTS. (`#159 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/159>`__)
- :func:`aiapy.calibrate.register` ``use_scipy`` keyword has been removed and replaced with a ``method`` keyword that defaults to use ``scipy`` by default.
  It is also possible to use ``scikit-image`` or ``opencv`` or ``cupy`` (provided you have either one installed.) (`#159 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/159>`__)


New Features
------------

- Added a "cupy" method to :func:`aiapy.calibrate.register` that will use cupy to do the affine_transform. (`#159 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/159>`__)


0.6.4 (2022-01-14)
==================

Bug Fixes
---------

- Fixes a bug where columns in the pointing table used to update the pointing information were being converted
  to masked float values.
  This bugfix ensures that any column used in the pointing update does not have a mask and any values that
  are masked are filled with NaN.
  This bug arises in astropy>=5.0. (`#151 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/151>`__)


Documentation
-------------

- Fixed escaping of characters in equations in multiple docstrings. (`#146 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/146>`__)


0.6.3 (2021-11-05)
==================

Bug Fixes
---------

- Fixes a bug in `aiapy.calibrate.update_pointing` concerning how the row in 3-hourly
  master pointing table is chosen.
  Previously, the row with ``T_START`` closest to ``DATE_OBS`` was chosen.
  Now, the row corresponding to ``T_OBS`` greater than or equal to ``T_START`` AND
  less than ``T_STOP`` is chosen. (`#137 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/137>`__)
- Update the ``x0_mp`` and ``y0_mp`` keywords when updating the pointing information
  in `aiapy.calibrate.update_pointing`. (`#140 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/140>`__)


Internal Changes
----------------

- In the case where a submap is passed into `aiapy.calibrate.fetch_spikes`,
  create the full-frame WCS directly from the submap WCS rather than creating
  an intermediate dummy full-frame map. (`#139 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/139>`__)


0.6.2 (2021-11-02)
==================

Bug Fixes
---------

- Fixed a bug in the units on the table returned by `aiapy.calibrate.util.get_pointing_table`.
  The ``X0`` and ``Y0`` columns were incorrectly being assigned units of arcseconds instead
  of pixels. (`#132 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/132>`__)
- Fixed an off-by-one bug in `aiapy.calibrate.update_pointing` where the
  ``CRPIX1`` and ``CRPIX2`` keywords were not being properly updated from the
  ``X0`` and ``Y0`` columns in the master pointing table. (`#132 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/132>`__)


0.6.1 (2021-11-01)
==================

Bug Fixes
---------

- Fixed a compatibility issue with sunpy>=3.1 in which creating a full-frame WCS in
  `aiapy.calibrate.fetch_spikes` was throwing an exception. (`#126 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/126>`__)
- Added a check on `aiapy.calibrate.update_pointing` so that passing in a submap or a map not at the
  full AIA resolution, raises an exception. (`#127 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/127>`__)


0.6.0 (2021-10-27)
==================

Breaking Changes
----------------

- Pin minimum version of Python to 3.7 (`#114 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/114>`__)
- Pin minimum version of sunpy to 3.0.0 LTS.

New Features
------------

- Added a new function :func:`aiapy.util.telescope_number` that returns the associated
  telscope number for a given filter wavelength. (`#116 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/116>`__)
- Added a new function :func:`aiapy.calibrate.util.get_error_table` to fetch and parse the
  table with the associate error parameters.
  This is used primarily in :func:`aiapy.calibrate.estimate_error`. (`#116 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/116>`__)
- Added a new function :func:`aiapy.calibrate.estimate_error` to calculate the error for
  a given set of AIA counts and the associated channel.
  This is an exact port of the IDL function ``aia_bp_estimate_error``. (`#116 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/116>`__)

Bug Fixes
---------

- `aiapy.calibrate.update_pointing` now skips updating keywords if the pointing values
  are missing from the pointing table returned from JSOC. (`#120 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/120>`__)

0.5.1 (2021-05-24)
==================

Backwards Incompatible Changes
------------------------------

- Pin sunpy dependency to ``<=3.0`` to allow for backwards compatibility with ``search_metadata``. (`#111 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/111>`__)

0.5.0 (2021-04-09)
==================

Features
--------

- Add a flag to :func:`aiapy.psf.deconvolve` that sets negative intensity values to zero before performing the deconvolution. (`#107 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/107>`__)

0.4.0 (2020-12-10)
==================

Features
--------

- Added a function (:func:`aiapy.util.sdo_location`) to obtain the SDO location at a given time. (`#57 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/57>`__)
- Added a function (:func:`aiapy.calibrate.respike`) for reinserting hot pixels into level 1 images. (`#62 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/62>`__)
- Updated default calibration version to 10.
  Added test for version 10 (`#90 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/90>`__)

Bug Fixes
---------

- Updated default calibration version number for degradation correction.
  Added tests for multiple calibration versions (`#74 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/74>`__)
- Fixed a bug where an out of date calibration epoch was used if there were older duplicate versions available in the same epoch. (`#90 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/90>`__)
- `aiapy.calibrate.util.get_pointing_table` now raises a more user-friendly `RuntimeError` if no pointing information can be found during the requested times.
  Previously it would raise a `KeyError`. (`#91 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/91>`__)
- `aiapy.calibrate.update_pointing` now searches 12 hours either side of the map date for pointing information.
  This allows for some very rare instances where more than 3 hours elapses between pointing information updates. (`#91 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/91>`__)

0.3.2 (2020-11-29)
==================

No significant changes.

0.3.1 (2020-11-15)
==================

Features
--------

- :func:`aiapy.calibrate.register` now raises a warning if the level number is missing or is greater than 1. (`#94 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/94>`__)

0.3.0 (2020-10-06)
==================

Features
--------

- Added a function (:func:`aiapy.calibrate.normalize_exposure`) to normalize an image by its exposure time. (`#78 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/78>`__)
- :func:`aiapy.calibrate.degradation` can now accept `~astropy.time.Time` objects with length greater than 1.
  This makes it easier to compute the channel degradation over long intervals. (`#80 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/80>`__)
- Citation information for `aiapy` is now available from `aiapy.__citation__`. (`#82 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/82>`__)
- The pointing table can now be passsed in as a keyword argument to :func:`aiapy.calibrate.update_pointing`.
  Added a :func:`aiapy.calibrate.util.get_pointing_table` to retrieve the 3-hour pointing table from JSOC over a given time interval. (`#84 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/84>`__)

Bug Fixes
---------

- The ``CROTA2`` keyword update in :func:`aiapy.calibrate.update_pointing` now includes the value of ``SAT_ROT`` from the FITS header.
  Previously, the keyword was only being updated with ``INSTROT``. (`#84 <https://gitlab.com/LMSAL_HUB/aia_hub/aiapy/-/merge_requests/84>`__)

0.2.0 (2020-07-16)
==================

Features
--------

- Functionality for respiking level 1 images and fetching spike data from JSOC
- Updated calibration data now fetched from JSOC to account for instrument degradation
- Compatibility fix with sunpy > 2.0.0 which previously caused level 1.5 maps to expand by several pixels
- Functionality for fetching the location of SDO in time

0.1.0  (2020-03-31)
===================

Features
--------

- Update pointing keywords in the header using the 3-hour pointing values from the JSOC
- Correct Heliographic Stonyhurst observer location
- Register images by removing the roll angle, centering the image, and scaling to a common resolution (i.e. "aia_prep")
- Calculate wavelength response functions for all channels, including time-dependent effects
- Account for channel degradation in image correction
- Compute the point spread function and deconvolve an image with the point spread function (with optional GPU acceleration)
