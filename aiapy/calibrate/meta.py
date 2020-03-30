"""
Functions for updating/fixing header keywords
"""
import copy

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import (SkyCoord, HeliocentricMeanEcliptic,
                                 CartesianRepresentation)
from sunpy.net import jsoc, attrs

__all__ = ['fix_observer_location', 'update_pointing']


def fix_observer_location(smap):
    """
    Fix inaccurate ``HGS_LON`` and ``HGS_LAT`` FITS keywords.

    The heliographic Stonyhurst latitude and longitude locations in the
    AIA FITS headers are incorrect. This function fixes the values of these
    keywords using the heliocentric aries ecliptic keywords, ``HAEX_OBS,
    HAEY_OBS, HAEZ_OBS``.

    .. note:: `~sunpy.map.sources.AIAMap` already accounts for the inaccurate
              HGS keywords by using the HAE keywords to construct the
              derived observer location.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
    """
    # Create observer coordinate from HAE coordinates
    coord = SkyCoord(
        x=smap.meta['haex_obs'] * u.m,
        y=smap.meta['haey_obs'] * u.m,
        z=smap.meta['haez_obs'] * u.m,
        representation_type=CartesianRepresentation,
        frame=HeliocentricMeanEcliptic,
        obstime=smap.date,
    ).heliographic_stonyhurst
    # Update header
    new_meta = copy.deepcopy(smap.meta)
    new_meta['hgln_obs'] = coord.lon.to(u.degree).value
    new_meta['hglt_obs'] = coord.lat.to(u.degree).value
    new_meta['dsun_obs'] = coord.radius.to(u.m).value

    return smap._new_instance(smap.data,
                              new_meta,
                              plot_settings=smap.plot_settings,
                              mask=smap.mask)


def update_pointing(smap):
    """
    Update map header to use the most recent 3-hourly pointing information
    from JSOC.

    This function queries JSOC for the 3-hour pointing table and updates the
    ``CRPIX1, CRPIX2, CDELT1, CDELT2, CROTA2`` keywords in the map header.

    Parameters
    ----------
    smap : `~sunpy.map.Map`

    Notes
    -----
    The method removes any ``PCi_j`` matrix keys in the header, and updates the
    ``CROTA2``.

    """
    # Query 3h pointing table from JSOC
    # NOTE: should this be a separate function?
    w_str = f'{smap.wavelength.to(u.angstrom).value:03.0f}'
    q = jsoc.JSOCClient().search_metadata(
        # Make range wide enough to get closest 3-hour pointing
        attrs.Time(smap.date - 3*u.h, end=smap.date + 3*u.h),
        attrs.jsoc.Series('aia.master_pointing3h'),
        attrs.jsoc.Keys(['T_START',
                         f'A_{w_str}_X0',
                         f'A_{w_str}_Y0',
                         f'A_{w_str}_INSTROT',
                         f'A_{w_str}_IMSCALE']),
    )
    # Find row closest to obstime
    table = Table.from_pandas(q)
    table['T_START'] = Time(table['T_START'], scale='utc')
    i_nearest = np.fabs((table['T_START'] - smap.date).to(u.s)).argmin()
    # Update headers
    new_meta = copy.deepcopy(smap.meta)
    new_meta['CRPIX1'] = table[f'A_{w_str}_X0'][i_nearest]
    new_meta['CRPIX2'] = table[f'A_{w_str}_Y0'][i_nearest]
    new_meta['CDELT1'] = table[f'A_{w_str}_IMSCALE'][i_nearest]
    new_meta['CDELT2'] = table[f'A_{w_str}_IMSCALE'][i_nearest]
    new_meta['CROTA2'] = table[f'A_{w_str}_INSTROT'][i_nearest]

    # SunPy map converts crota to a PCi_j matrix, so we remove it to force the
    # re-conversion.
    new_meta.pop('PC1_1')
    new_meta.pop('PC1_2')
    new_meta.pop('PC2_1')
    new_meta.pop('PC2_2')

    return smap._new_instance(smap.data,
                              new_meta,
                              plot_settings=smap.plot_settings,
                              mask=smap.mask)
