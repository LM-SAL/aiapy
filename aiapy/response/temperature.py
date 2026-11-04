"""
Temperature response functions for the AIA EUV channels.
"""

import numpy as np
from sunkit_instruments.response import SourceSpectra

import astropy.units as u
from astropy.constants import c, h
from astropy.io import fits

from sunpy.time import parse_time

import aiapy.calibrate
from aiapy.data._manager import manager

__all__ = ["get_source_spectra", "temperature_response"]

_EMISSIVITY_URL = "https://raw.githubusercontent.com/LM-SAL/backup-files/refs/heads/aia/static/aiapy/"
_EMISSIVITY_FILES = {"solarsoft": "emissivity_v9", "chianti11": "emissivity_chianti11"}


@manager.require(
    "emissivity_v9",
    [_EMISSIVITY_URL + "aia_V9_emissivity.fits"],
    "b926fe19f7daabd596f5b42d8997114223ebf7a5f434ba09d086a9f8cf2bd24b",
    defer_download=True,
)
@manager.require(
    "emissivity_chianti11",
    [_EMISSIVITY_URL + "aia_chianti11_emissivity.fits"],
    "90dbd44f455a8a6e9c2f0aca2f32474ad5dd8ea8167851aad6f751ed546c252f",
    defer_download=True,
)
def get_source_spectra(source="solarsoft"):
    r"""
    The emissivity of an optically thin plasma for the AIA temperature responses.

    Both sources are summed line and continuum emissivities on a 0.1 Å grid from 1 to 901 Å
    and log T from 4 to 9 in steps of 0.05, with the plasma of SolarSoft's
    ``aia_V9_fullemiss.genx``: Feldman 1992 extended coronal abundances, the CHIANTI
    ionization equilibrium, a constant pressure of :math:`10^{15}` K cm\\ :sup:`-3`,
    and the He II 304 Å doublet boosted by 20 for its optical depth.
    SolarSoft's V9 file boosts only the weak member of the doublet, because its exact
    wavelength match for the strong one failed on CHIANTI 9 (V6 and V8 boosted both);
    ``"chianti11"`` boosts both, so its 304 Å response is about 2.7 times SolarSoft's
    at He II temperatures.
    ``meta`` records these settings, and for ``"solarsoft"`` the ``aia_V9_chiantifix.genx``
    correction used by `temperature_response`.

    Parameters
    ----------
    source : {"solarsoft", "chianti11"}, optional
        ``"solarsoft"`` is the emissivity that SolarSoft uses, from ``aia_V9_fullemiss.genx``
        (CHIANTI 9.0.1), with each line moved to the nearest grid wavelength.
        SolarSoft converts it to DN with :math:`hc` = 12398 eV Å, so it is rescaled here by
        12398 eV Å / :math:`hc`; folding it with the exact :math:`hc` that
        `aiapy.response.Channel` uses then reproduces SolarSoft.
        ``"chianti11"`` is computed from CHIANTI 11.0.2 with fiasco, with each line shared
        between the two grid wavelengths around it, so folding it with a linearly interpolated
        effective area gives the response at the line's exact wavelength.
        Between CHIANTI 9 and 11 the Fe XIV line at 93.614 Å moved to 91.273 Å, out of the
        94 Å passband, which lowers the 94 Å response at log T 6.3 by two thirds; the
        chiantifix does not apply to it.

    Returns
    -------
    `sunkit_instruments.response.SourceSpectra`
    """
    with fits.open(manager.get(_EMISSIVITY_FILES[source])) as hdul:
        header = hdul[0].header
        # astype gives native byte order, which xarray needs
        temperature = 10 ** hdul["TEMPERATURE"].data["logte"].astype(float) * u.K
        wavelength = hdul["WAVELENGTH"].data["wavelength"].astype(float) * u.angstrom
        emissivity = hdul["EMISSIVITY"].data.astype(float) * u.Unit("ph cm3 s-1 sr-1 Angstrom-1")
        meta = {
            "chianti_version": header["CHIANTI"],
            "abundance": header["ABUNDANC"],
            "ioneq": header["IONEQ"],
            "pressure": header["PRESSURE"] * u.K / u.cm**3,
            "fixes": list(header["HISTORY"]),
        }
        if "EMISFILE" in header:  # SolarSoft's file: its hc and its chiantifix
            emissivity *= (12398 * u.eV * u.angstrom / (h * c)).to_value(u.dimensionless_unscaled)
            chiantifix = hdul["CHIANTIFIX"].data
            meta["source"] = header["EMISFILE"]
            meta["chiantifix"] = {
                str(channel): row.astype(float) * u.Unit("cm5 DN / (pix s)")
                for channel, row in zip(chiantifix["channel"], chiantifix["chiantifix"], strict=True)
            }
        else:
            meta["source"] = f"fiasco {header['FIASCO']}"
    return SourceSpectra(temperature, wavelength, emissivity, density=meta["pressure"] / temperature, meta=meta)


def temperature_response(channel, spectra, obstime=None, *, chiantifix=False):
    r"""
    Temperature response of an AIA EUV channel.

    The response is

    .. math::

        R(T, t) = \left[R_0(T) C_e + R_{fix}(T)\right] C_T(t),

    where :math:`R_0 C_e` folds ``spectra`` with the wavelength response of ``channel``
    (:math:`C_e` is its EVE correction, 1 unless ``include_eve_correction=True``),
    :math:`R_{fix}` is the SolarSoft chiantifix and :math:`C_T` is the time-dependent
    degradation from `aiapy.calibrate.degradation`.
    With ``get_source_spectra()`` this is SolarSoft's
    ``aia_get_response(/temperature, /dn)``, with ``/evenorm`` if the channel has
    ``include_eve_correction=True`` and ``timedepend_date`` if ``obstime`` is given.
    The results agree with SolarSoft to about :math:`10^{-4}`.

    The chiantifix is SolarSoft's empirical correction for 94 Å emission missing from
    CHIANTI; it is zero for the other channels.
    It raises the warm 94 Å response by up to a factor of 2.7 near log T = 6.2 and
    lowers it by about a quarter above log T = 6.4.
    It was derived for CHIANTI 7.1.3 and the EVE-normalized response and has not been
    updated since, so it is only applied to the SolarSoft emissivity with
    ``include_eve_correction=True``.
    Its effect on six-channel DEMs is small, but analyses that use 94 Å on its own
    depend on it.

    Parameters
    ----------
    channel : `aiapy.response.Channel`
        An EUV channel.
    spectra : `sunkit_instruments.response.SourceSpectra`
        Source spectra, for example from `get_source_spectra`.
    obstime : `~astropy.time.Time`, optional
        Time or times of the observation.
        If `None` (the default), no time-dependent degradation is applied.
    chiantifix : `bool`, optional
        If true, add the SolarSoft chiantifix.

    Returns
    -------
    temperature : `~astropy.units.Quantity`
    response : `~astropy.units.Quantity`
        In cm5 DN / (pix s). An array ``obstime`` adds a leading time axis.
    """
    response = spectra.temperature_response(channel)
    if chiantifix:
        if "chiantifix" not in spectra.meta:
            msg = 'chiantifix only applies to the SolarSoft emissivity from get_source_spectra("solarsoft")'
            raise ValueError(msg)
        if not channel.include_eve_correction:
            msg = "chiantifix was derived for the EVE-normalized response; use include_eve_correction=True"
            raise ValueError(msg)
        response = response + spectra.meta["chiantifix"][channel.name]
    if obstime is not None:
        obstime = parse_time(obstime)
        degradation = aiapy.calibrate.degradation(
            channel.channel,
            obstime,
            correction_table=channel._correction_table,
            calibration_version=channel._calibration_version,
        )
        response = degradation[:, np.newaxis] * response
        if obstime.isscalar:
            response = response[0]
    return spectra.temperature, response
