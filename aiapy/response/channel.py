"""
Class for accessing response function data from each channel
"""
import os
import collections

import numpy as np
import astropy.units as u
import astropy.constants as const
import astropy.time
import astropy.table
from sunpy.io.special import read_genx
from sunpy.util.config import get_and_create_download_dir
from sunpy.util.net import check_download_file
from sunpy.util.metadata import MetaDict

from .util import get_correction_table

__all__ = ['Channel']

SSW_AIA_URL = 'https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/response/'
AIA_INSTRUMENT_FILE = 'aia_V{}_{}_fullinst.genx'  # What changes with version?
VERSION_NUMBER = 8  # Most recent version number for instrument data


class Channel(object):
    """
    Object for accessing AIA channel properties and calculating response
    functions

    Parameters
    ----------
    channel : `~astropy.units.Quantity`
        Wavelength of AIA channel
    instrument_file : `str`, optional
        Path to AIA instrument file. If not specified, the latest version will
        be downloaded from SolarSoft.
    """

    @u.quantity_input
    def __init__(self, channel: u.angstrom, instrument_file=None):
        self._channel = channel
        self._instrument_data = self._get_instrument_data(instrument_file)

    @property
    def is_fuv(self):
        """
        Returns True for UV and visible channels: 1600, 1700, 4500 |AA|.

        .. |AA| unicode:: x212B .. angstrom
        """
        return self.channel in [1600, 1700, 4500]*u.angstrom

    def _get_instrument_data(self, instrument_file):
        """
        Read the raw instrument data for all channels from the `.genx` files
        in SSW
        """
        # TODO: Do not rely on genx files for this data!
        # Avoid reading data if possible
        if isinstance(instrument_file, collections.OrderedDict):
            return instrument_file
        # TODO: refactor using the new remote data capabilities in sunpy
        if instrument_file is None:
            # TODO: don't put this in $HOME/sunpy
            local_dir = get_and_create_download_dir()
            filename = AIA_INSTRUMENT_FILE.format(
                VERSION_NUMBER, 'fuv' if self.is_fuv else 'all')
            check_download_file(filename, SSW_AIA_URL, local_dir)
            instrument_file = os.path.join(local_dir, filename)
        return read_genx(instrument_file)

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
    def plate_scale(self,) -> u.pixel / u.steradian:
        return u.Quantity(self._data['platescale'], u.pixel/u.steradian)

    @property
    @u.quantity_input
    def effective_area(self,) -> u.cm**2:
        """
        Uncorrected effective area as a function of wavelength

        According to Section 2 of [1]_, the effective area is given by,

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
        .. [1] Boerner et al., 2012, Sol. Phys., `275, 41 <http://adsabs.harvard.edu/abs/2012SoPh..275...41B>`_
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

    def _select_epoch_from_table(self, obstime, **kwargs):
        """
        Return correction table with only the first epoch and the epoch in
        which `obstime` falls.
        """
        if kwargs.get('correction_table', None) is not None:
            table = kwargs.get('correction_table')
        else:
            table = get_correction_table()
        # Select only this channel
        table = table[table['WAVE_STR'] == f'{self.name}_THIN']
        # Select only most recent version number, JSOC keeps some old entries
        table = table[table['VER_NUM'] == VERSION_NUMBER]
        # Select the epoch for the given observation time
        obstime_in_epoch = np.logical_and(obstime >= table['T_START'],
                                          obstime < table['T_STOP'])
        if not obstime_in_epoch.any():
            raise IndexError(f'No valid calibration epoch for {obstime}')
        # Create new table with only first and obstime epochs
        return astropy.table.Table(rows=[table[0], table[obstime_in_epoch][0]],
                                   names=table.colnames)

    @u.quantity_input
    def time_correction(self, obstime, **kwargs) -> u.dimensionless_unscaled:
        """
        Correction to effective area to account for time-dependent degradation
        of instrument

        The correction factor to account for the time-varying degradation of
        the telescopes is given by a normalization to the calibration epoch
        closest to `obstime` and an interpolation within that epoch to
        `obstime`,

        .. math::

            \\frac{A_{eff}(t_{e})}{A_{eff}(t_0)}(1 + p_1\delta t + p_2\delta t^2 + p_3\delta t^3)

        where :math:`A_{eff}(t_e)` is the effective area calculated at the
        calibration epoch for `obstime`, :math:`A_{eff}(t_0)` is the effective
        area at the first calibration epoch (i.e. at launch),
        :math:`p_1,p_2,p_3` are the interpolation coefficients for the
        `obstime` epoch, and :math:`\delta t` is the difference between the
        start time of the epoch and `obstime`.

        All calibration terms are taken from the `aia.response` series in JSOC
        or read from the table input by the user. This function is adapted
        directly from the `aia_bp_corrections.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_corrections.pro>`_ routine in SolarSoft.

        Parameters
        ----------
        obstime : `~astropy.time.Time`

        See Also
        --------
        wavelength_response
        eve_correction
        """
        table = self._select_epoch_from_table(obstime, **kwargs)
        # Time difference between obstime and start of epoch
        dt = (obstime - table['T_START'][-1]).to(u.day).value
        # Correction to most recent epoch
        ratio = table['EFF_AREA'][-1] / table['EFF_AREA'][0]
        # Polynomial correction to interpolate within epoch
        poly = (table['EFFA_P1'][-1]*dt
                + table['EFFA_P2'][-1]*dt**2
                + table['EFFA_P3'][-1]*dt**3
                + 1.)
        return u.Quantity(poly * ratio)

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
        the `aia_bp_corrections.pro <https://hesperia.gsfc.nasa.gov/ssw/sdo/aia/idl/response/aia_bp_corrections.pro>`_ routine in SolarSoft.

        Parameters
        ----------
        obstime : `~astropy.time.Time`
        """
        table = self._select_epoch_from_table(obstime, **kwargs)
        effective_area_interp = np.interp(table['EFF_WVLN'][-1],
                                          self.wavelength.to(u.angstrom),
                                          self.effective_area.to(u.cm**2))
        return u.Quantity(table['EFF_AREA'][0] / effective_area_interp)

    @property
    @u.quantity_input
    def gain(self,) -> u.count / u.ph:
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
        correction_table : `~astropy.table.Table`, optional
            Table of correction parameters. If not specified, it will be
            queried from JSOC.

        See Also
        --------
        effective_area : Uncorrected effective area as a function of wavelength
        gain : Gain of the telescope for the particular channel
        crosstalk : Correction factor for crosstalk channel
        eve_correction : Correction factor from EVE calibration
        time_correction : Correction factor for time-dependent degradation
        """
        eve_correction, time_correction = 1, 1
        if obstime is not None:
            time_correction = self.time_correction(obstime, **kwargs)
            if include_eve_correction:
                eve_correction = self.eve_correction(obstime, **kwargs)
        crosstalk = self.crosstalk if include_crosstalk else 0*u.cm**2
        return ((self.effective_area + crosstalk)
                * self.gain
                * time_correction
                * eve_correction)

    @u.quantity_input
    def temperature_response(self,) -> u.count / u.pixel / u.s * u.cm**5:
        raise NotImplementedError(
            'Temperature response calculation has not yet been implemented')
