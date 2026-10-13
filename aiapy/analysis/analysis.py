"""
Functions for analyzing AIA images.
"""

import numpy as np

import astropy.units as u

# Warren et al. (2012), Eq. A1, lists four coefficients but inconsistent bounds;
# use a cubic with coefficients in ascending power order.
_FE18_COEFFICIENTS = (-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3)
# 171 Å weight in the 171/193 Å composite.
_FE18_MIXTURE_WEIGHT = 0.31
# Median composite and 94 Å intensities of the fit data, in DN/s.
_FE18_COMPOSITE_SCALE = 116.54
_FE18_INTENSITY_SCALE = 0.39
# The fit poorly constrains scaled composite intensities above 30.
_FE18_COMPOSITE_LIMIT = 30.0


def isolate_fe18_intensity(map_1, map_2, map_3):
    r"""
    Estimate Fe XVIII 93.9 Å emission in AIA 94 Å observations.

    Subtract the warm 94 Å contribution estimated from the 171 Å and 193 Å
    channels using the empirical fit of :cite:t:`warren2012`.

    .. warning::

        These fixed coefficients were fitted to 2010 March 22 observations.
        They may not suit other dates or calibrations; check the warm estimate
        for your data.

    Parameters
    ----------
    map_1, map_2, map_3 : `~sunpy.map.GenericMap`
        Level 1.5 (degraded corrected optional but recommended) images at
        94, 171, and 193 Å, in any order and in DN/s.
        They must share a pixel grid and be close enough in time for the source.
        Units, wavelengths, and dimensions are checked; alignment and times are not.

    Returns
    -------
    `~sunpy.map.GenericMap`
        The 94 Å image minus the estimated warm contribution, with all input
        masks combined.

    Notes
    -----
    The warm contribution is estimated as

    .. math::

        I_{94,\mathrm{warm}} = 0.39 \, \mathrm{DN/s}\sum_{i=0}^{3} a_i
        \left[\frac{f I_{171} + (1 - f) I_{193}}{116.54 \, \mathrm{DN/s}}\right]^i

    where :math:`a = [-7.31\times10^{-2}, 9.75\times10^{-1}, 9.90\times10^{-2},
    -2.84\times10^{-3}]`, :math:`f = 0.31`, and 0.39 and 116.54 DN/s are the
    median intensities of the 94 Å and composite fit images. The scaled
    composite intensity is capped at 30 because brighter fit data are sparse.
    Normalize each input by its exposure time, e.g. ``smap / smap.exposure_time``.

    Equation (A1) has inconsistent bounds; its four coefficients define powers
    0 to 3 here.

    This method does not apply to flares, where Fe XXIV 192.04 Å enters the
    193 Å channel, and is unreliable for bright 1 MK emission such as moss.
    """
    maps = (map_1, map_2, map_3)
    for name, smap in zip(("map_1", "map_2", "map_3"), maps, strict=True):
        if smap.unit != u.DN / u.s:
            msg = (
                f"{name} must have units of DN/s, got {smap.unit}. "
                "Level 1.5 images should be divided by their exposure times, "
                "e.g. smap / smap.exposure_time."
            )
            raise ValueError(msg)
    wavelengths = [smap.wavelength.to_value(u.angstrom) if smap.wavelength is not None else None for smap in maps]
    if set(wavelengths) != {94, 171, 193}:
        msg = f"Expected one map each at 94, 171, and 193 Å; got {wavelengths} Å"
        raise ValueError(msg)
    m_94, m_171, m_193 = (smap for _, smap in sorted(zip(wavelengths, maps, strict=True)))
    if not (m_94.data.shape == m_171.data.shape == m_193.data.shape):
        msg = (
            "All images must have the same dimensions, got "
            f"{m_94.data.shape}, {m_171.data.shape}, and {m_193.data.shape}"
        )
        raise ValueError(msg)
    composite = (_FE18_MIXTURE_WEIGHT * m_171.data + (1 - _FE18_MIXTURE_WEIGHT) * m_193.data) / _FE18_COMPOSITE_SCALE
    composite = np.minimum(composite, _FE18_COMPOSITE_LIMIT)
    warm = _FE18_INTENSITY_SCALE * np.polynomial.polynomial.polyval(composite, _FE18_COEFFICIENTS)
    result = m_94 - warm * m_94.unit
    masks = [np.broadcast_to(smap.mask, m_94.data.shape) for smap in (m_94, m_171, m_193) if smap.mask is not None]
    if masks:
        result.mask = np.logical_or.reduce(masks)
    return result
