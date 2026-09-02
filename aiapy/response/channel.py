"""
Class for accessing response function data from each channel.
"""

from pathlib import Path
from urllib.parse import urljoin
from collections.abc import Mapping

import numpy as np
from sunkit_instruments.response.abstractions import AbstractChannel

import astropy.units as u

from sunpy.io.special import read_genx
from sunpy.util.metadata import MetaDict

import aiapy.calibrate
from aiapy import _SSW_MIRRORS
from aiapy.calibrate.utils import LATEST_CORRECTION_VERSION, _select_epoch_from_correction_table, get_correction_table
from aiapy.data._manager import manager
from aiapy.utils import telescope_number
from aiapy.utils.decorators import validate_channel

__all__ = ["Channel"]

AIA_INSTRUMENT_FILE = "sdo/aia/response/aia_V{}_{}_fullinst.genx"
# Most recent version number for instrument response data, there is 9 but its the same as V8.
VERSION_NUMBER = 8
URL_HASH = {
    8: {
        "fuv": (
            [urljoin(mirror, AIA_INSTRUMENT_FILE.format(VERSION_NUMBER, "fuv")) for mirror in _SSW_MIRRORS],
            "8635166d8f6dde48da4f135925f4e8f48a0574f129c2c2ca24da6628550f5430",
        ),
        "euv": (
            [urljoin(mirror, AIA_INSTRUMENT_FILE.format(VERSION_NUMBER, "all")) for mirror in _SSW_MIRRORS],
            "3940648e6b02876c45a9893f40806bbcc50baa994ae3fa2d95148916988426dd",
        ),
    },
}


class Channel(AbstractChannel):
    """
    Interface to AIA channel properties and response functions.

    This class provides an interface to the AIA channels and methods
    for calculating the effective area and wavelength response functions
    as a function of wavelength.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of AIA channel.
    instrument_file : `str` or `dict`, optional
        Path to AIA instrument file or the already-parsed contents of one
        (as returned by `~sunpy.io.special.read_genx`).
        If not specified, the latest version will be downloaded from SolarSoft.
    include_eve_correction : `bool`, optional
        If true, include correction to EVE calibration in effective area.
    include_crosstalk : `bool`, optional
        If true, include the effect of crosstalk between channels that share a telescope
    correction_table : `~astropy.table.Table` or `str`, optional
        Table of correction parameters or path to correction table file.
        If not specified, the latest version is fetched from SSW.
    calibration_version : `int`, optional
        The version of the calibration to use when calculating the
        degradation. By default, this is the most recent version available
        from JSOC. If you are using a specific calibration response file,
        you may need to specify this according to the version in that file.

    Examples
    --------
    >>> import astropy.units as u
    >>> from aiapy.response import Channel
    >>> c = Channel(171 * u.angstrom)  # doctest: +REMOTE_DATA
    >>> c.telescope_number  # doctest: +REMOTE_DATA
    3
    >>> c.name  # doctest: +REMOTE_DATA
    '171'
    >>> c.channel  # doctest: +REMOTE_DATA
    <Quantity 171. Angstrom>
    """

    @u.quantity_input
    @validate_channel("channel")
    def __init__(
        self,
        channel: u.angstrom,
        *,
        instrument_file=None,
        include_eve_correction=False,
        include_crosstalk=True,
        correction_table=None,
        calibration_version=None,
    ) -> None:
        self._channel = channel
        self._instrument_data = self._get_instrument_data(instrument_file)
        self.include_eve_correction = include_eve_correction
        self.include_crosstalk = include_crosstalk
        self.correction_table = correction_table
        self.calibration_version = calibration_version

    @property
    def is_fuv(self):
        """
        Returns True for UV and visible channels 1600, 1700, 4500 Å.
        """
        return self.channel in [1600, 1700, 4500] * u.angstrom

    def _get_instrument_data(self, instrument_file):
        """
        Read the raw instrument data for all channels from the ``.genx`` files
        in SSW.
        """
        if isinstance(instrument_file, Mapping):
            return instrument_file
        if instrument_file is None:
            instrument_file = self._get_fuv_instrument_file() if self.is_fuv else self._get_euv_instrument_file()
        return read_genx(instrument_file)

    @manager.require("instrument_file_euv", *URL_HASH[VERSION_NUMBER]["euv"])
    def _get_euv_instrument_file(self):
        return manager.get("instrument_file_euv")

    @manager.require("instrument_file_fuv", *URL_HASH[VERSION_NUMBER]["fuv"])
    def _get_fuv_instrument_file(self):
        return manager.get("instrument_file_fuv")

    @property
    def _data(
        self,
    ):
        """
        Instrument data for this channel.
        """
        return MetaDict(self._instrument_data[f"A{self.name}_FULL"])

    @property
    @u.quantity_input
    def channel(
        self,
    ) -> u.angstrom:
        """
        Nominal wavelength at which the bandpass of the channel is centered.
        """
        return self._channel

    @property
    def name(
        self,
    ) -> str:
        return f"{self.channel.to(u.angstrom).value:.0f}"

    @property
    def telescope_number(
        self,
    ):
        """
        Label denoting the telescope to which the given channel is assigned.

        See `crosstalk` for context of why this is important.
        """
        return telescope_number(self.channel)

    @property
    @u.quantity_input
    def wavelength(
        self,
    ) -> u.angstrom:
        """
        Array of wavelengths over which channel properties are calculated.
        """
        return u.Quantity(self._data["wave"], u.angstrom)

    @property
    @u.quantity_input
    def primary_mirror_reflectance(
        self,
    ) -> u.dimensionless_unscaled:
        return u.Quantity(self._data["primary"])

    @property
    @u.quantity_input
    def secondary_mirror_reflectance(
        self,
    ) -> u.dimensionless_unscaled:
        return u.Quantity(self._data["secondary"])

    @property
    @u.quantity_input
    def mirror_reflectance(self):
        """
        Combined reflectance of the primary and secondary mirrors.
        """
        return self.primary_mirror_reflectance * self.secondary_mirror_reflectance

    @property
    @u.quantity_input
    def focal_plane_filter_transmittance(
        self,
    ) -> u.dimensionless_unscaled:
        return u.Quantity(self._data["fp_filter"])

    @property
    @u.quantity_input
    def entrance_filter_transmittance(
        self,
    ) -> u.dimensionless_unscaled:
        return u.Quantity(self._data["ent_filter"])

    @property
    @u.quantity_input
    def filter_transmittance(self):
        """Combined transmittance of the focal plane and entrance filters"""
        return self.entrance_filter_transmittance * self.focal_plane_filter_transmittance

    @property
    @u.quantity_input
    def geometrical_area(
        self,
    ) -> u.cm**2:
        return u.Quantity(self._data["geoarea"], u.cm**2)

    @property
    @u.quantity_input
    def effective_quantum_efficiency(
        self,
    ) -> u.dimensionless_unscaled:
        return u.Quantity(self._data["ccd"])

    @property
    @u.quantity_input
    def _preflight_contamination(self):
        # NOTE: Preflight contamination calibration data missing for FUV channels
        return u.Quantity(self._data.get("contam", np.ones(self.wavelength.shape)))

    @property
    @u.quantity_input
    def energy_per_electron(self) -> u.eV / u.electron:
        return 1 / self._data["elecperev"] * u.eV / u.electron

    @property
    @u.quantity_input
    def camera_gain(self) -> u.DN / u.electron:
        return 1 / self._data["elecperdn"] * u.DN / u.electron

    @property
    @u.quantity_input
    def pixel_solid_angle(
        self,
    ) -> u.steradian / u.pixel:
        return u.Quantity(self._data["platescale"], u.steradian / u.pixel)

    @property
    def correction_table(self):
        """
        Correction table as provided: a table, a path, or `None`.
        """
        return self._correction_table_input

    @correction_table.setter
    def correction_table(self, value):
        self._correction_table_input = value
        self._resolved_correction_table = None

    @property
    def _correction_table(self):
        """
        The correction table resolved to a table: read/fetched when the
        ``correction_table`` attribute is a path or `None` and cached until
        ``correction_table`` is set again.
        """
        if self._resolved_correction_table is None:
            source = self.correction_table
            if source is None:
                self._resolved_correction_table = get_correction_table()
            elif isinstance(source, str | Path):
                self._resolved_correction_table = get_correction_table(source=source)
            else:
                self._resolved_correction_table = source
        return self._resolved_correction_table

    @property
    def _calibration_version(self):
        return LATEST_CORRECTION_VERSION if self.calibration_version is None else self.calibration_version

    @u.quantity_input
    def _get_eve_correction(self, obstime) -> u.dimensionless_unscaled:
        """
        Calculate degradation factor from EVE rocket flight cross-calibration.

        This function is adapted directly from the
        `aia_bp_corrections.pro <https://sohoftp.nascom.nasa.gov/solarsoft/sdo/aia/idl/response/aia_bp_corrections.pro>`__
        routine in SolarSoft.
        """
        table = _select_epoch_from_correction_table(
            self.channel,
            obstime,
            self._correction_table,
            calibration_version=self._calibration_version,
        )
        # NOTE: Explicitly interpolating to the pristine (time-independent,
        # crosstalk-free) effective area since this correction is part of the
        # degradation. This matches aia_bp_corrections.pro, which uses the
        # effective area as stored in the instrument file.
        effective_area_interp = np.interp(table["EFF_WVLN"][-1], self.wavelength, self._pristine_effective_area)
        return table["EFF_AREA"][0] / effective_area_interp

    @u.quantity_input
    def degradation(self, obstime=None) -> u.dimensionless_unscaled:
        r"""
        Time-dependent degradation factor.

        The AIA time-dependent degradation is modeled as,

        .. math::

            D(t) = C_e(t)C_d(t)

        where :math:`C_e(t)` is the EVE cross-calibration correction factor
        and :math:`C_d(t)` is the time-varying telescope degradation computed
        by `aiapy.calibrate.degradation`. The time-independent pre-flight
        contamination model :math:`D_c(\lambda)` (Section 2.1.6 of [boerner])
        is part of the static effective area, see `effective_area`.

        The EVE correction factor is given by,

        .. math::

            C_e(t) = \frac{A_{eff}(\lambda_n,t_0)}{A_{eff}(\lambda_E,t_e)}

        where :math:`A_{eff}(\lambda_n,t_0)` is the effective area at the
        nominal wavelength of the channel (:math:`\lambda_n`) at the first
        calibration epoch and :math:`A_{eff}(\lambda_E,t_e)` is the effective
        area at the ```obstime``` calibration epoch interpolated to the effective
        wavelength (:math:`\lambda_E`).

        Parameters
        ----------
        obstime : `~astropy.time.Time`, optional
            The time of the observation. If `None` (the default), no
            degradation is applied.

        Returns
        -------
        `~astropy.units.Quantity`

        See Also
        --------
        aiapy.calibrate.degradation

        References
        ----------
        .. [boerner] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`__
        """
        if obstime is None:
            return u.Quantity(1.0)
        degradation = aiapy.calibrate.degradation(
            self.channel,
            obstime,
            correction_table=self._correction_table,
            calibration_version=self._calibration_version,
        )
        if self.include_eve_correction:
            degradation = degradation * self._get_eve_correction(obstime)
        return degradation

    @property
    @u.quantity_input
    def _pristine_effective_area(self) -> u.cm**2:
        """
        Nominal effective area with no crosstalk or time-dependent
        degradation.
        """
        return (
            self.geometrical_area
            * self.mirror_reflectance
            * self.filter_transmittance
            * self.effective_quantum_efficiency
            * self._preflight_contamination
        )

    @property
    @u.quantity_input
    def _crosstalk(self) -> u.cm**2:
        """
        Pristine cross-talk component of the effective area.
        """
        crosstalk_lookup = {
            94 * u.angstrom: 304 * u.angstrom,
            304 * u.angstrom: 94 * u.angstrom,
            131 * u.angstrom: 335 * u.angstrom,
            335 * u.angstrom: 131 * u.angstrom,
        }
        if self.channel not in crosstalk_lookup:
            return u.Quantity(np.zeros(self.wavelength.shape), u.cm**2)
        cross = type(self)(crosstalk_lookup[self.channel], instrument_file=self._instrument_data)
        return (
            cross.geometrical_area
            * cross.mirror_reflectance
            * self.focal_plane_filter_transmittance
            * cross.entrance_filter_transmittance
            * cross.effective_quantum_efficiency
            * cross._preflight_contamination
        )

    @u.quantity_input
    def effective_area(self) -> u.cm**2:
        r"""
        Effective area as a function of wavelength.

        According to Section 2 of [boerner]_, the effective area is given by,

        .. math::

            A_{eff}(\lambda) = A_{geo}R_P(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)Q(\lambda)D_c(\lambda)D(t),

        where,

        - :math:`A_{geo}`: geometrical collecting area
        - :math:`R_P`, :math:`R_S`: reflectance of primary and secondary mirrors, respectively
        - :math:`T_E`, :math:`T_F`: transmittance of the entrance and focal-plane filters, respectively
        - :math:`Q`: effective quantum efficiency of the CCD
        - :math:`D_c`: preflight contaminant transmittance of the optics
        - :math:`D`: time-dependent degradation, evaluated at the time bound
          with `~sunkit_instruments.response.abstractions.AbstractChannel.at`
          (no time-dependent degradation for an unbound channel)

        The effective area contains information about the efficiency of the telescope optics and its sensitivity
        as a function of wavelength.
        All of the telescope properties are read from the AIA instrument files available in SolarSoft.

        On telescopes 1, 3, and 4, both channels are always illuminated.
        This can lead to contamination in a channel from the channel with which it shares a telescope.
        This impacts the 94 and 304 Å channels as well as 131 and 335 Å.
        For these channels, an additional component is added to the effective area to model the cross-talk
        between these channels.
        See Section 2.2.1 of [boerner]_ for more details.

        See Also
        --------
        degradation

        References
        ----------
        .. [boerner] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`__
        """
        effective_area = self._pristine_effective_area
        if self.include_crosstalk:
            effective_area = effective_area + self._crosstalk
        # NOTE: The nominal channel's time-dependent degradation is applied to
        # the total (nominal + crosstalk) effective area since the postflight
        # cross-calibration is done on the whole channel including crosstalk.
        if self._obstime is not None:
            effective_area = effective_area * self.degradation(obstime=self._obstime)
        return effective_area
