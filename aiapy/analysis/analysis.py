"""
Functions for analyzing AIA images.
"""

import numpy as np

import astropy.units as u

import sunpy.map

# Coefficients of the polynomial fit in Equation (A1) of the appendix in
# Warren et al. (2012), in order of increasing power of the composite intensity.
# Equation (A1) has inconsistent summation bounds; its four coefficients
# are interpreted here as a cubic, with powers from 0 to 3.
_FE18_COEFFICIENTS = (-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3)
# Weight of the 171 A image in the composite 171/193 A intensity.
_FE18_MIXTURE_WEIGHT = 0.31
# Scaling factors derived from the median intensities of the fit data,
# in units of DN/s.
_FE18_COMPOSITE_SCALE = 116.54
_FE18_INTENSITY_SCALE = 0.39
# Limit on the composite intensity, in units of the scaled composite intensity.
# Higher intensities are poorly constrained by the fit data.
_FE18_COMPOSITE_LIMIT = 30.0


def isolate_fe18_intensity(m_94, m_171, m_193):
    r"""
    Isolate the Fe XVIII emission in AIA 94 Å observations.

    The 94 Å channel is contaminated by "warm" emission formed at temperatures
    well below the roughly 7 MK at which Fe XVIII is formed. This function
    estimates that warm contribution from an empirical polynomial fit to a
    mixture of the 171 Å and 193 Å channels and subtracts it from the 94 Å
    channel, following the appendix of :cite:t:`warren2012`.

    Parameters
    ----------
    m_94 : `~sunpy.map.GenericMap`
        Co-aligned level 1.5 image of the 94 Å channel in units of DN/s.
    m_171 : `~sunpy.map.GenericMap`
        Co-aligned level 1.5 image of the 171 Å channel in units of DN/s.
    m_193 : `~sunpy.map.GenericMap`
        Co-aligned level 1.5 image of the 193 Å channel in units of DN/s.

    Returns
    -------
    `~sunpy.map.GenericMap`
        The 94 Å image with the estimated warm contribution subtracted,
        providing an estimate of the Fe XVIII 93.9 Å line intensity.
        Pixels masked in any input are masked in the output.

    Notes
    -----
    The warm contribution is estimated as

    .. math::

        I_{94,\mathrm{warm}} = 0.39 \, \mathrm{DN/s}\sum_{i=0}^{3} a_i
        \left[\frac{f I_{171} + (1 - f) I_{193}}{116.54 \, \mathrm{DN/s}}\right]^i

    where :math:`a = [-7.31\times10^{-2}, 9.75\times10^{-1}, 9.90\times10^{-2},
    -2.84\times10^{-3}]`, :math:`f = 0.31`, and 0.39 and 116.54 DN/s are the
    median intensities of the 94 Å channel and of the composite 171/193 Å
    image used for the fit. The scaled composite intensity is limited to 30,
    beyond which the fit is poorly constrained by the data. The intensities of
    the input images must be normalized by their exposure times, e.g. by
    dividing them by ``smap.exposure_time``.

    Equation (A1) has inconsistent summation bounds in the published and
    arXiv versions. Here, the four listed coefficients define a cubic with
    powers from 0 to 3.

    This procedure is not applicable to flare observations, where Fe XXIV
    192.04 Å contributes to the 193 Å channel, or to regions of very bright
    1 MK emission, such as moss.

    """
    for name, smap in (("m_94", m_94), ("m_171", m_171), ("m_193", m_193)):
        if not isinstance(smap, sunpy.map.GenericMap):
            msg = f"{name} must be a sunpy Map, got {type(smap)}"
            raise TypeError(msg)
        if smap.unit != u.DN / u.s:
            msg = (
                f"{name} must have units of DN/s, got {smap.unit}. "
                "Level 1.5 images should be divided by their exposure times, "
                "e.g. smap / smap.exposure_time."
            )
            raise ValueError(msg)
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
