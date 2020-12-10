Aiapy 0.4.0 (2020-12-10)
========================

Features
--------

- Added a function (:func:`~aiapy.util.sdo_location`) to obtain the SDO location at a given time. (`#57 <https://github.com/sunpy/aiapy/pull/57>`__)
- Added a function (:func:`~aiapy.calibrate.respike`) for reinserting hot pixels into level 1 images. (`#62 <https://github.com/sunpy/aiapy/pull/62>`__)
- Added a function (:func:`~aiapy.calibrate.normalize_exposure`) to normalize an image
  by its exposure time. (`#78 <https://github.com/sunpy/aiapy/pull/78>`__)
- :func:`~aiapy.calibrate.degradation` can now accept `~astropy.time.Time` objects with
  length greater than 1. This makes it easier to compute the channel degradation over
  long intervals. (`#80 <https://github.com/sunpy/aiapy/pull/80>`__)
- Citation information for `aiapy` is now available from `~aiapy.__citation__`. (`#82 <https://github.com/sunpy/aiapy/pull/82>`__)
- The pointing table can now be passsed in as a keyword argument to :func:`~aiapy.calibrate.update_pointing`.
  Added a :func:`~aiapy.calibrate.util.get_pointing_table` to retrieve the 3-hour pointing table
  from JSOC over a given time interval. (`#84 <https://github.com/sunpy/aiapy/pull/84>`__)
- Updated default calibration version to 10
  Added test for version 10 (`#90 <https://github.com/sunpy/aiapy/pull/90>`__)
- `aiapy.calibrate.register` now raises a warning if the level number
  is missing or is greater than 1. (`#94 <https://github.com/sunpy/aiapy/pull/94>`__)


Bug Fixes
---------

- Updated default calibration version number for degradation correction
  Added tests for multiple calibration versions (`#74 <https://github.com/sunpy/aiapy/pull/74>`__)
- The `CROTA2` keyword update in :func:`~aiapy.calibrate.update_pointing` now includes
  the value of `SAT_ROT` from the FITS header. Previously, the keyword was only being
  updated with `INSTROT`. (`#84 <https://github.com/sunpy/aiapy/pull/84>`__)
- Fixed a bug where an out of date calibration epoch was used if there were older
  duplicate versions available in the same epoch. (`#90 <https://github.com/sunpy/aiapy/pull/90>`__)
- `aiapy.calibrate.util.get_pointing_table` now raises a more user-friendly
  `RuntimeError` if no pointing information can be found during the requested
  times. Previously it would raise a `KeyError`. (`#91 <https://github.com/sunpy/aiapy/pull/91>`__)
- `aiapy.calibrate.update_pointing` now searches 12 hours either side of the map
  date for pointing information. This allows for some very rare instances where
  more than 3 hours elapses between pointing information updates. (`#91 <https://github.com/sunpy/aiapy/pull/91>`__)
