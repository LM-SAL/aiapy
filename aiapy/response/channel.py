"""
Class for accessing response function data from each channel
"""
import collections

import numpy as np
import astropy.units as u
import astropy.constants as const
from sunpy.io.special import read_genx
from sunpy.util.metadata import MetaDict
from sunpy.data import manager

from aiapy.calibrate.util import _select_epoch_from_table, get_correction_table
from aiapy.calibrate import degradation

__all__ = ['Channel']

AIA_INSTRUMENT_FILE = 'https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/aia_V{}_{}_fullinst.genx'  # What changes with version?
VERSION_NUMBER = 8  # Most recent version number for instrument response data
# URLs and SHA-256 hashes for each version for the EUV and FUV files
# The URLs are left as a list so that possible mirrors for these files
# can be specified
URL_HASH = {
    2: {'fuv': None, 'euv': None},
    3: {'fuv': None, 'euv': None},
    4: {'fuv': None, 'euv': None},
    6: {'fuv': None, 'euv': None},
    8: {'fuv': ((AIA_INSTRUMENT_FILE.format(VERSION_NUMBER, 'fuv')),
                '8635166d8f6dde48da4f135925f4e8f48a0574f129c2c2ca24da6628550f5430'),
        'euv': ((AIA_INSTRUMENT_FILE.format(VERSION_NUMBER, 'all'),),
                '3940648e6b02876c45a9893f40806bbcc50baa994ae3fa2d95148916988426dd')},
}


class Channel(object):
    """
    Interface to AIA channel properties and response functions.

    This class provides an interface to the AIA channels and methods
    for calculating the effective area and wavelength response functions
    as a function of wavelength.

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of AIA channel
    instrument_file : `str`, optional
        Path to AIA instrument file. If not specified, the latest version will
        be downloaded from SolarSoft.

    Examples
    ---------
    >>> import astropy.units as u
    >>> from aiapy.response import Channel
    >>> c = Channel(171*u.angstrom)  # doctest: +REMOTE_DATA
    >>> c.telescope_number  # doctest: +REMOTE_DATA
    3
    >>> c.name  # doctest: +REMOTE_DATA
    '171'
    >>> c.channel  # doctest: +REMOTE_DATA
    <Quantity 171. Angstrom>
    """

    @u.quantity_input
    def __init__(self, channel: u.angstrom, instrument_file=None):
        self._channel = channel
        self._instrument_data = self._get_instrument_data(instrument_file)

    @property
    def is_fuv(self):
        """
        Returns True for UV and visible channels 1600, 1700, 4500 |AA|.

        .. |AA| unicode:: x212B .. angstrom
        """
        return self.channel in [1600, 1700, 4500]*u.angstrom

    def _get_instrument_data(self, instrument_file):
        """
        Read the raw instrument data for all channels from the `.genx` files
        in SSW
        """
        if isinstance(instrument_file, collections.OrderedDict):
            return instrument_file
        if instrument_file is None:
            if self.is_fuv:
                instrument_file = self._get_fuv_instrument_file()
            else:
                instrument_file = self._get_euv_instrument_file()
        return read_genx(instrument_file)

    @manager.require('instrument_file_euv', *URL_HASH[VERSION_NUMBER]['euv'])
    def _get_euv_instrument_file(self):
        return manager.get('instrument_file_euv')

    @manager.require('instrument_file_fuv', *URL_HASH[VERSION_NUMBER]['fuv'])
    def _get_fuv_instrument_file(self):
        return manager.get('instrument_file_fuv')

    @property
    def _data(self,):
        """
        Instrument data for this channel
        """
        return MetaDict(self._instrument_data[f'A{self.name}_FULL'])

    @property
    @u.quantity_input
    def channel(self,) -> u.angstrom:
        """
        Nominal wavelength at which the bandpass of the channel is centered
        """
        return self._channel

    @property
    def name(self,):
        return f'{self.channel.to(u.angstrom).value:.0f}'

    @property
    def telescope_number(self,):
        """
        Label denoting the telescope to which the given channel is assigned.
        See `crosstalk` for context of why this is important.
        """
        return {
            94*u.angstrom: 4,
            131*u.angstrom: 1,
            171*u.angstrom: 3,
            193*u.angstrom: 2,
            211*u.angstrom: 2,
            304*u.angstrom: 4,
            335*u.angstrom: 1,
            1600*u.angstrom: 3,
            1700*u.angstrom: 3,
            4500*u.angstrom: 3,
        }[self.channel]

    @property
    @u.quantity_input
    def wavelength(self,) -> u.angstrom:
        """
        Array of wavelengths over which channel properties are calculated
        """
        return u.Quantity(self._data['wave'], u.angstrom)

    @property
    @u.quantity_input
    def primary_reflectance(self,) -> u.dimensionless_unscaled:
        return u.Quantity(self._data['primary'])

    @property
    @u.quantity_input
    def secondary_reflectance(self,) -> u.dimensionless_unscaled:
        return u.Quantity(self._data['secondary'])

    @property
    @u.quantity_input
    def focal_plane_filter_efficiency(self,) -> u.dimensionless_unscaled:
        return u.Quantity(self._data['fp_filter'])

    @property
    @u.quantity_input
    def entrance_filter_efficiency(self,) -> u.dimensionless_unscaled:
        return u.Quantity(self._data['ent_filter'])

    @property
    @u.quantity_input
    def geometrical_collecting_area(self,) -> u.cm**2:
        return u.Quantity(self._data['geoarea'], u.cm**2)

    @property
    @u.quantity_input
    def quantum_efficiency(self,) -> u.dimensionless_unscaled:
        return u.Quantity(self._data['ccd'])

    @property
    @u.quantity_input
    def contamination(self,) -> u.dimensionless_unscaled:
        # Contamination missing for FUV channels
        if 'contam' in self._data:
            return u.Quantity(self._data['contam'])
        else:
            return u.Quantity([1])

    @property
    @u.quantity_input
    def plate_scale(self,) -> u.steradian / u.pixel:
        return u.Quantity(self._data['platescale'], u.steradian/u.pixel)

    @property
    @u.quantity_input
    def effective_area(self,) -> u.cm**2:
        """
        Uncorrected effective area as a function of wavelength

        According to Section 2 of [boerner]_, the effective area is given by,

        .. math::

            A_{eff}(\lambda) = A_{geo}R_P(\lambda)R_S(\lambda)T_E(\lambda)T_F(\lambda)D(\lambda)Q(\lambda),

        where,

        - :math:`A_{geo}`: geometrical collecting area
        - :math:`R_P`, :math:`R_S`: reflectances of primary and secondary
          mirrors, respectively
        - :math:`T_E`, :math:`T_F`: transmission efficiency of the entrance
          and focal-plane filters, respectively
        - :math:`D`: contaminant transmittance of optics
        - :math:`Q`: quantum efficiency of the CCD

        The effective area contains information about the efficiency of the
        telescope optics and its sensitivity as a function of wavelength. All
        of the telescope properties are read from the AIA instrument files
        available in SolarSoft.

        References
        ----------
        .. [boerner] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`_
        """
        return (self.primary_reflectance
                * self.secondary_reflectance
                * self.focal_plane_filter_efficiency
                * self.entrance_filter_efficiency
                * self.geometrical_collecting_area
                * self.quantum_efficiency
                * self.contamination)

    @property
    @u.quantity_input
    def crosstalk(self,) -> u.cm**2:
        """
        Contamination of effective area from crosstalk  between channels.

        On telescopes 1, 3, and 4, both channels are always illuminated. This
        can lead to contamination in a channel from the channel with which it
        shares a telescope. This impacts the 94 and 304 channels as well as
        131 and 335. See Section 2.2.1 of [1]_ for more details.

        References
        ----------
        .. [1] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`_
        """
        crosstalk_lookup = {
            94*u.angstrom: 304*u.angstrom,
            304*u.angstrom: 94*u.angstrom,
            131*u.angstrom: 335*u.angstrom,
            335*u.angstrom: 131*u.angstrom,
        }
        if self.channel in crosstalk_lookup:
            cross = Channel(crosstalk_lookup[self.channel],
                            instrument_file=self._instrument_data)
            return (cross.primary_reflectance
                    * cross.secondary_reflectance
                    * self.focal_plane_filter_efficiency
                    * cross.entrance_filter_efficiency
                    * cross.geometrical_collecting_area
                    * cross.quantum_efficiency
                    * cross.contamination)
        else:
            return u.Quantity(np.zeros(self.wavelength.shape), u.cm**2)

    @u.quantity_input
    def eve_correction(self, obstime, **kwargs) -> u.dimensionless_unscaled:
        """
        Correct effective area to give good agreement with full-disk EVE data

        The EVE correction factor is given by,

        .. math::

            \\frac{A_{eff}(\lambda_n,t_0)}{A_{eff}(\lambda_E,t_e)}

        where :math:`A_{eff}(\lambda_n,t_0)` is the effective area at the
        nominal wavelength of the channel (:math:`\lambda_n`) at the first
        calibration epoch and :math:`A_{eff}(\lambda_E,t_e)` is the effective
        area at the `obstime` calibration epoch interpolated to the effective
        wavelength (:math:`\lambda_E`). This function is adapted directly from
        the `aia_bp_corrections.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_corrections.pro>`_
        routine in SolarSoft.

        Parameters
        ----------
        obstime : `~astropy.time.Time`

        Other Parameters
        ---------------------
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

        Returns
        -------
        `~astropy.units.Quantity`

        See Also
        --------
        aiapy.calibrate.util.get_correction_table
        """
        table = _select_epoch_from_table(
            self.channel,
            obstime,
            get_correction_table(correction_table=kwargs.get('correction_table')),
            version=kwargs.get('calibration_version'),
        )
        effective_area_interp = np.interp(table['EFF_WVLN'][-1],
                                          self.wavelength,
                                          self.effective_area)
        return table['EFF_AREA'][0] / effective_area_interp

    @property
    @u.quantity_input
    def gain(self,) -> u.count / u.ph:
        """
        Gain of the CCD camera system.

        According to Section 2 of [boerner1]_, the gain of the CCD-camera system, in
        DN per photon, is given by,

        .. math::

            G(\lambda) = \\frac{hc}{\lambda}\\frac{g}{a}

        where :math:`g` is the camera gain in DN per electron
        and :math:`a\\approx 3.65` eV per electron is a conversion factor.

        References
        ----------
        .. [boerner1] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`_
        """
        _e = u.electron  # Avoid rewriting u.electron a lot
        electron_per_ev = self._data['elecperev'] * _e / u.eV
        energy_per_photon = const.h * const.c / self.wavelength / u.ph
        electron_per_photon = (electron_per_ev
                               * energy_per_photon).to(_e / u.ph)
        # Cannot discharge less than one electron per photon
        discharge_floor = electron_per_photon < (1 * _e / u.ph)
        electron_per_photon[discharge_floor] = 1 * _e / u.ph
        return electron_per_photon / (self._data['elecperdn'] * _e / u.count)

    @u.quantity_input
    def wavelength_response(self,
                            obstime=None,
                            include_eve_correction=False,
                            include_crosstalk=True,
                            **kwargs) -> u.count / u.ph * u.cm**2:
        """
        The wavelength response function is the product of the gain and the
        effective area.

        The wavelength response as a function of time and wavelength is given
        by,

        .. math::

            R(\lambda,t) = (A_{eff}(\lambda) + A_{cross}(\lambda))G(\lambda)C_T(t)C_E(t)

        where,

        - :math:`A_{eff}(\lambda)` is the effective area as a function of
          wavelength
        - :math:`A_{cross}(\lambda)` is the effective area of the crosstalk
          channel
        - :math:`G(\lambda)` is the gain of the telescope
        - :math:`C_T(t)` is the time-dependent correction factor for the
          instrument degradation
        - :math:`C_E(t)` is the time-dependent EVE correction factor

        Parameters
        ----------
        obstime : `~astropy.time.Time`, optional
            If specified, a time-dependent correction is applied to account
            for degradation
        include_eve_correction : `bool`, optional
            If true and `obstime` is not `None`, include correction to EVE
            calibration. The time-dependent correction is also included.
        include_crosstalk : `bool`, optional
            If true, include the effect of crosstalk between channels that
            share a telescope

        Other Parameters
        ----------------
        correction_table : `~astropy.table.Table` or `str`, optional
            Table of correction parameters or path to correction table file.
            If not specified, it will be queried from JSOC. See
            `~aiapy.calibrate.util.get_correction_table` for more information.

        Returns
        -------
        `~astropy.units.Quantity`

        See Also
        --------
        effective_area
        gain
        crosstalk
        eve_correction
        aiapy.calibrate.degradation
        """
        eve_correction, time_correction = 1, 1
        if obstime is not None:
            time_correction = degradation(self.channel, obstime, **kwargs)
            if include_eve_correction:
                eve_correction = self.eve_correction(obstime, **kwargs)
        crosstalk = self.crosstalk if include_crosstalk else 0*u.cm**2
        return ((self.effective_area + crosstalk)
                * self.gain
                * time_correction
                * eve_correction)
