Aiapy 0.3.1 (2020-11-15)
========================

Features
--------

- `aiapy.calibrate.register` now raises a warning if the level number
  is missing or is greater than 1. (`#94 <https://github.com/sunpy/aiapy/pull/94>`__)


aiapy v0.3.0 (2020-10-06)
=========================

Features
--------

- Added a function (:func:`~aiapy.calibrate.normalize_exposure`) to normalize an image
  by its exposure time. (`#78 <https://github.com/sunpy/aiapy/pull/78>`__)
- :func:`~aiapy.calibrate.degradation` can now accept `~astropy.time.Time` objects with
  length greater than 1. This makes it easier to compute the channel degradation over
  long intervals. (`#80 <https://github.com/sunpy/aiapy/pull/80>`__)
- Citation information for `aiapy` is now available from `~aiapy.__citation__`. (`#82 <https://github.com/sunpy/aiapy/pull/82>`__)
- The pointing table can now be passsed in as a keyword argument to :func:`~aiapy.calibrate.update_pointing`.
  Added a :func:`~aiapy.calibrate.util.get_pointing_table` to retrieve the 3-hour pointing table
  from JSOC over a given time interval. (`#84 <https://github.com/sunpy/aiapy/pull/84>`__)


Bug Fixes
---------

- The `CROTA2` keyword update in :func:`~aiapy.calibrate.update_pointing` now includes
  the value of `SAT_ROT` from the FITS header. Previously, the keyword was only being
  updated with `INSTROT`. (`#84 <https://github.com/sunpy/aiapy/pull/84>`__)
