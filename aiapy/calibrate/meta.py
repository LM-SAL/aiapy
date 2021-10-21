"""
Functions for updating/fixing header keywords
"""
import copy
import warnings

import numpy as np

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic, SkyCoord

from aiapy.calibrate.util import get_pointing_table
from aiapy.util.exceptions import AiapyUserWarning

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
    smap : `~sunpy.map.source.sdo.AIAMap`
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

    return smap._new_instance(
        smap.data,
        new_meta,
        plot_settings=smap.plot_settings,
        mask=smap.mask
    )


def update_pointing(smap, pointing_table=None):
    """
    Update pointing information in the map header

    This function updates the pointing information in the map by
    updating the ``CRPIX1, CRPIX2, CDELT1, CDELT2, CROTA2`` keywords
    in the header using the information provided in `pointing_table`.
    If `pointing_table` is not specified, the 3-hour pointing
    information is queried from `JSOC <http://jsoc.stanford.edu/>`_.

    .. note:: The method removes any ``PCi_j`` matrix keys in the header and
              updates the ``CROTA2`` keyword.

    .. note:: If correcting pointing information for a large number of images,
              it is strongly recommended to query the table once for the
              appropriate interval and then pass this table in rather than
              executing repeated queries.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
    pointing_table : `~astropy.table.QTable`, optional
        Table of pointing information. If not specified, the table
        will be retrieved from JSOC.

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap`

    See Also
    --------
    aiapy.calibrate.util.get_pointing_table
    """
    if pointing_table is None:
        # Make range wide enough to get closest 3-hour pointing
        pointing_table = get_pointing_table(smap.date - 12*u.h, smap.date + 12*u.h)
    # Find row closest to obstime
    i_nearest = np.fabs((pointing_table['T_START'] - smap.date).to(u.s)).argmin()
    w_str = f'{smap.wavelength.to(u.angstrom).value:03.0f}'
    new_meta = copy.deepcopy(smap.meta)
    # Extract new pointing parameters
    crpix1 = pointing_table[f'A_{w_str}_X0'][i_nearest].to('arcsecond').value
    crpix2 = pointing_table[f'A_{w_str}_Y0'][i_nearest].to('arcsecond').value
    cdelt = pointing_table[f'A_{w_str}_IMSCALE'][i_nearest].to('arcsecond / pixel').value
    # CROTA2 is the sum of INSTROT and SAT_ROT.
    # See http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_H_AIA-SDO_FITS_Keyword_Document.pdf
    # NOTE: Is the value of SAT_ROT in the header accurate?
    crota2 = pointing_table[f'A_{w_str}_INSTROT'][i_nearest] + smap.meta['SAT_ROT'] * u.degree
    crota2 = crota2.to('deg').value
    # Update headers
    for key, value in [('crpix1', crpix1),
                       ('crpix2', crpix2),
                       ('cdelt1', cdelt),
                       ('cdelt2', cdelt),
                       ('crota2', crota2)]:
        if np.isnan(value):
            # There are some entries in the pointing table returned from the JSOC that are marked as
            # MISSING. These get converted to NaNs when we cast it to an astropy quantity table. In
            # these cases, we just want to skip updating the pointing information.
            warnings.warn(f'Missing value in pointing table for {key}. This key will not be updated.',
                          AiapyUserWarning)
        else:
            new_meta[key] = value

    # sunpy map converts crota to a PCi_j matrix, so we remove it to force the
    # re-conversion.
    new_meta.pop('PC1_1')
    new_meta.pop('PC1_2')
    new_meta.pop('PC2_1')
    new_meta.pop('PC2_2')
    return smap._new_instance(
        smap.data,
        new_meta,
        plot_settings=smap.plot_settings,
        mask=smap.mask
    )
