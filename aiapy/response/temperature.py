"""
Temperature response functions for the AIA EUV channels.
"""

import numpy as np
from sunkit_instruments.response import SourceSpectra

import astropy.units as u
from astropy.io import fits

from sunpy.time import parse_time

import aiapy.calibrate
from aiapy.data._manager import manager

__all__ = ["get_source_spectra", "temperature_response"]

_EMISSIVITY_URL = (
    "https://raw.githubusercontent.com/LM-SAL/backup-files/refs/heads/aia/static/aiapy/aia_V9_emissivity.fits"
)
_EMISSIVITY_SHA256 = "b926fe19f7daabd596f5b42d8997114223ebf7a5f434ba09d086a9f8cf2bd24b"


@manager.require("emissivity_v9", [_EMISSIVITY_URL], _EMISSIVITY_SHA256)
def get_source_spectra():
    r"""
    The emissivity that SolarSoft uses for the AIA temperature responses.

    This is the summed line and continuum emissivity of ``aia_V9_fullemiss.genx``
    (CHIANTI 9.0.1, Feldman 1992 extended coronal abundances, ``chianti.ioneq``,
    constant pressure of :math:`10^{15}` K cm\\ :sup:`-3`, He II 303.786 Å boosted by 20),
    on a 0.1 Å grid from 1 to 901 Å.
    ``meta`` records these settings and the ``aia_V9_chiantifix.genx`` correction
    used by `temperature_response`.

    Returns
    -------
    `sunkit_instruments.response.SourceSpectra`
    """
    with fits.open(manager.get("emissivity_v9")) as hdul:
        header = hdul[0].header
        # astype gives native byte order, which xarray needs
        temperature = 10 ** hdul["TEMPERATURE"].data["logte"].astype(float) * u.K
        wavelength = hdul["WAVELENGTH"].data["wavelength"].astype(float) * u.angstrom
        emissivity = hdul["EMISSIVITY"].data.astype(float) * u.Unit("ph cm3 s-1 sr-1 Angstrom-1")
        chiantifix = hdul["CHIANTIFIX"].data
        meta = {
            "chianti_version": header["CHIANTI"],
            "abundance": header["ABUNDANC"],
            "ioneq": header["IONEQ"],
            "pressure": header["PRESSURE"] * u.K / u.cm**3,
            "fixes": list(header["HISTORY"]),
            "source": header["EMISFILE"],
            "chiantifix": {
                str(channel): row.astype(float) * u.Unit("cm5 DN / (pix s)")
                for channel, row in zip(chiantifix["channel"], chiantifix["chiantifix"], strict=True)
            },
        }
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
    The results agree with SolarSoft to about :math:`10^{-4}`; :math:`3.4\times10^{-5}`
    of that is because SolarSoft takes :math:`hc` as 12398 eV Å.

    The chiantifix is SolarSoft's empirical correction for 94 Å emission missing from
    CHIANTI; it is zero for the other channels.
    It was derived for CHIANTI 7.1.3 and the EVE-normalized response and has not been
    updated since, so it is only applied to the SolarSoft emissivity with
    ``include_eve_correction=True``. It changes the 94 Å response by up to a factor of
    2.7, at log T = 6.2.

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
            msg = "chiantifix only applies to the SolarSoft emissivity from get_source_spectra()"
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
