"""
==================================================
Correcting automatic exposure control light curves
==================================================

During flares, AIA interleaves long fixed exposures (AEC off) with short
exposures set by the automatic exposure control (AEC on).
Once divided by the exposure time, both should give the same count rate, but
the AEC-on light curve shows a sawtooth that follows the exposure-time ramp.
This can look like quasi-periodic pulsations that are not real.

.. warning::

    The coefficients used here were tuned by hand for this flare.
    They may not suit other observations; compare the AEC-on and AEC-off light curves for your data.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import astropy.units as u
from astropy.time import Time

from sunpy.net import Fido
from sunpy.net import attrs as a

from aiapy.calibrate import correct_aec

plt.rcParams["date.converter"] = "concise"

###############################################################################
# The correction only needs header keywords, so we query the JSOC for the
# metadata of the 131 Å images and do not download any data.

results = Fido.search(
    a.Time("2014-09-08T23:00:00", "2014-09-09T04:00:00"),
    a.jsoc.Series.aia_lev1_euv_12s,
    a.Wavelength(131 * u.angstrom),
)
table = results["jsoc"].show("T_OBS", "AECTYPE", "DATAMEAN", "EXPTIME", "NSATPIX")
print(table)

###############################################################################
# Without any coefficients, `~aiapy.calibrate.correct_aec` returns
# ``DATAMEAN / EXPTIME``.
# The corrected count rate accounts for a dark pedestal that is 0.084 DN too
# low, AEC-on exposures that are 3.2 ms longer than ``EXPTIME`` and residual
# charge of 1 DN/s for every 18000 saturated pixels in the neighbouring images.

raw = correct_aec(table)
corrected = correct_aec(
    table,
    dn_bias=-0.084 * u.DN,
    exptime_bias=-3.2 * u.ms,
    saturation_scale=18000,
    saturation_shift=0 * u.s,
)

###############################################################################
# To compare the two AEC states, we interpolate the AEC-off count rate to the
# times of the AEC-on images and take the difference.

time = Time(table["T_OBS"])
seconds = (time - time[0]).to_value(u.s)
time = time.datetime64
aec_on = table["AECTYPE"] > 1


def aec_difference(rate):
    """
    AEC-on count rate minus the AEC-off count rate interpolated to the same times.
    """
    return rate[aec_on] - np.interp(seconds[aec_on], seconds[~aec_on], rate[~aec_on])


raw_difference = aec_difference(raw)
corrected_difference = aec_difference(corrected)
print(
    f"Mean absolute difference: uncorrected {np.mean(np.abs(raw_difference)):.4f}, corrected {np.mean(np.abs(corrected_difference)):.4f}"
)

###############################################################################
# The top panel shows the uncorrected light curves, where the AEC-on one has a
# sawtooth that follows its exposure time.
# The bottom panel shows that the sawtooth is mostly gone after the correction.

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, height_ratios=[3, 2], layout="constrained")
ax1.plot(time[~aec_on], raw[~aec_on], label="AEC off", color="#52514e")
ax1.plot(time[aec_on], raw[aec_on], label="AEC on", color="#2a78d6")
ax1.set_ylabel(f"Mean count rate [{raw.unit:latex}]")
ax1.legend()
ax2.axhline(0, color="#d0d0d0", linewidth=1)
ax2.plot(time[aec_on], raw_difference, label="uncorrected", color="#2a78d6")
ax2.plot(time[aec_on], corrected_difference, label="corrected", color="#eb6834")
ax2.set_ylabel(f"AEC on - AEC off [{raw.unit:latex}]")
ax2.set_xlabel("Time (UTC)")
ax2.legend()

plt.show()

###############################################################################
# For other observations, the coefficients can be fitted by minimizing the
# difference between the two AEC states.
# The pedestal and exposure-time errors both grow as ``EXPTIME`` gets shorter,
# so a fit struggles to separate them. Instead, measure the pedestal and fit
# only the other two.


def mismatch(params):
    """
    Mean absolute difference between the corrected AEC-on and AEC-off count rates.
    """
    rate = correct_aec(
        table,
        dn_bias=-0.084 * u.DN,
        exptime_bias=params[0] * u.ms,
        saturation_scale=1e4 * params[1],
        saturation_shift=0 * u.s,
    )
    return np.mean(np.abs(aec_difference(rate))).value


fit = minimize(mismatch, x0=[0, 1], method="Nelder-Mead")
print(
    f"exptime_bias: {fit.x[0]:.1f} ms, saturation_scale: {1e4 * fit.x[1]:.0f}, mean absolute difference: {fit.fun:.4f}"
)

###############################################################################
# The same works for your own light curve, such as the mean of a region.
# Build the table from the level 1 files of every frame, with the mean of the
# region in DN (not divided by the exposure time) as ``DATAMEAN``, and measure
# the pedestal error from a dark corner of the images:
#
# .. code-block:: python
#
#     maps = sunpy.map.Map(files)
#     table = Table(
#         {
#             "T_OBS": [m.meta["t_obs"] for m in maps],
#             "AECTYPE": [m.meta["aectype"] for m in maps],
#             "DATAMEAN": [m.submap(bottom_left, top_right=top_right).data.mean() for m in maps],
#             "EXPTIME": [m.meta["exptime"] for m in maps],
#             "NSATPIX": [m.meta["nsatpix"] for m in maps],
#         }
#     )
#     dn_bias = np.mean([m.data[:64, :64].mean() for m in maps]) * u.DN
#
# ``NSATPIX`` counts the saturated pixels in the whole image, while
# ``saturation_scale`` is relative to the mean of your region, so fit it again
# for each region.
