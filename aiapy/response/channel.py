"""
Class for accessing response function data from each channel.
"""

import collections
from urllib.parse import urljoin

import astropy.units as u
import numpy as np
from sunkit_instruments.response.abstractions import AbstractChannel
from sunpy.data import manager
from sunpy.io.special import read_genx
from sunpy.util.metadata import MetaDict

import aiapy.calibrate
from aiapy import _SSW_MIRRORS
from aiapy.calibrate.util import _select_epoch_from_correction_table, get_correction_table
from aiapy.util import telescope_number
from aiapy.util.decorators import validate_channel

__all__ = ["Channel"]

# TODO: Work out what changes with version.
AIA_INSTRUMENT_FILE = "sdo/aia/response/aia_V{}_{}_fullinst.genx"
VERSION_NUMBER = 8  # Most recent version number for instrument response data
# URLs and SHA-256 hashes for each version for the EUV and FUV files
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
    instrument_file : `str`, optional
        Path to AIA instrument file.
        If not specified, the latest version will be downloaded from SolarSoft.
    include_eve_correction : `bool`, optional
        If true, include correction to EVE calibration in effective area.
    include_crosstalk : `bool`, optional
        If true, include the effect of crosstalk between channels that share a telescope
    correction_table : `~astropy.table.Table` or `str`, optional
        Table of correction parameters or path to correction table file.
        If not specified, it will be queried from JSOC.
        If you are calling this function repeatedly, it is recommended to
        read the correction table once and pass it with this argument to avoid
        multiple redundant network calls.
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
        **kwargs,
    ):
        self._channel = channel
        self._instrument_data = self._get_instrument_data(instrument_file)
        self.include_eve_correction = include_eve_correction
        self.include_crosstalk = include_crosstalk
        self.correction_table = kwargs.get("correction_table", None)
        self.calibration_version = kwargs.get("calibration_version", None)

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
        if isinstance(instrument_file, collections.OrderedDict):
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
    ):
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
        "Combined transmittance of the focal plane and entrance filters"
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
            get_correction_table(correction_table=self.correction_table),
            version=self.calibration_version,
        )
        # NOTE: Explicitly setting obstime=None here because we should be
        # interpolating to the uncorrected EA since this is part of the
        # correction.
        effective_area_interp = np.interp(table["EFF_WVLN"][-1], self.wavelength, self.effective_area(obstime=None))
        return table["EFF_AREA"][0] / effective_area_interp

    @u.quantity_input
    def degradation(self, obstime=None) -> u.dimensionless_unscaled:
        r"""
        Wavelength- and time-dependent degradation factor.

        The AIA wavelength- and time-dependent degradation is modeled as,

        .. math::

            D(\lambda, t) = D_c(\lambda)C_e(t)C_d(t)

        where :math:`D_c(\lambda)` is the pre-flight contamination model
        as described in Section 2.1.6 of [boerner],
        :math:`C_e(t)` is the EVE cross-calibration correction factor,
        and :math:`C_d(t)` is the time-varying telescope degradation computed
        by `aiapy.calibrate.degradation`.

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
        obstime : `~astropy.time.Time`
            The time of the observation.

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
        contamination = self._preflight_contamination
        if obstime is not None:
            contamination *= aiapy.calibrate.degradation(
                self.channel,
                obstime,
                correction_table=self.correction_table,
                calibration_version=self.calibration_version,
            )
            if self.include_eve_correction:
                contamination *= self._get_eve_correction(obstime)
        return contamination

    @u.quantity_input
    def _get_crosstalk(self, obstime) -> u.cm**2:
        """
        Compute cross-talk component of effective area.
        """
        crosstalk_lookup = {
            94 * u.angstrom: 304 * u.angstrom,
            304 * u.angstrom: 94 * u.angstrom,
            131 * u.angstrom: 335 * u.angstrom,
            335 * u.angstrom: 131 * u.angstrom,
        }
        if self.channel not in crosstalk_lookup:
            return u.Quantity(np.zeros(self.wavelength.shape), u.cm**2)
        cross = type(self)(
            crosstalk_lookup[self.channel],
            instrument_file=self._instrument_data,
            include_eve_correction=self.include_eve_correction,
            include_crosstalk=self.include_crosstalk,
            correction_table=self.correction_table,
            calibration_version=self.calibration_version,
        )
        effective_area = (
            cross.geometrical_area
            * cross.mirror_reflectance
            * self.focal_plane_filter_transmittance
            * cross.entrance_filter_transmittance
            * cross.effective_quantum_efficiency
            * cross._preflight_contamination  # noqa: SLF001
        )
        # NOTE: This applies only the time-dependent component of the degradation for
        # the nominal channel. The logic here is that the nominal time-dependent correction
        # should be applied to the total (nominal + crosstalk) effective area since the
        # cross-calibration is done on the whole channel postflight which includes the crosstalk
        effective_area *= self.degradation(obstime=obstime) / self._preflight_contamination
        return effective_area

    @u.quantity_input
    def effective_area(
        self,
        *,
        obstime=None,
    ) -> u.cm**2:
        r"""
        Effective area as a function of wavelength.

        According to Section 2 of [boerner]_, the effective area is given by,

        .. math::

            A_{eff}(\lambda) = A_{geo}R_P(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)Q(\lambda)D(\lambda,t),

        where,

        - :math:`A_{geo}`: geometrical collecting area
        - :math:`R_P`, :math:`R_S`: reflectance of primary and secondary mirrors, respectively
        - :math:`T_E`, :math:`T_F`: transmissitance of the entrance and focal-plane filters, respectively
        - :math:`Q`: effective quantum efficiency of the CCD
        - :math:`D`: degradation of optics, including preflight contaminant transmittance and time-dependent corrections

        The effective area contains information about the efficiency of the telescope optics and its sensitivity
        as a function of wavelength.
        All of the telescope properties are read from the AIA instrument files available in SolarSoft.

        On telescopes 1, 3, and 4, both channels are always illuminated.
        This can lead to contamination in a channel from the channel with which it shares a telescope.
        This impacts the 94 and 304 Å channels as well as 131 and 335 Å.
        For these channels, additional component is added to the effective area to model the cross-talk
        between these channels.
        See Section 2.2.1 of [boerner]_ for more details.

        See Also
        --------
        degradation

        References
        ----------
        .. [boerner] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`__
        """
        effective_area = super().effective_area(obstime=obstime)
        if self.include_crosstalk:
            effective_area += self._get_crosstalk(obstime)
        return effective_area
