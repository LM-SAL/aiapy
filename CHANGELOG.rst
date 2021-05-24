Aiapy 0.5.1 (2021-05-24)
========================

Backwards Incompatible Changes
------------------------------

- Pin sunpy dependency to ``<=3.0`` to allow for backwards compatibility with ``search_metadata``. (`#111 <https://github.com/sunpy/aiapy/pull/111>`__)


Aiapy 0.5.0 (2021-04-09)
========================

Features
--------

- Add a flag to :func:`aiapy.psf.deconvolve` that sets negative intensity values
  to zero before performing the deconvolution. (`#107 <https://github.com/sunpy/aiapy/pull/107>`__)
