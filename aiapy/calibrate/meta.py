"""
Functions for updating/fixing header keywords.
"""
import copy
import warnings

import numpy as np

import astropy.time
import astropy.units as u
from astropy.coordinates import CartesianRepresentation, HeliocentricMeanEcliptic, SkyCoord
from sunpy.map import contains_full_disk

from aiapy.calibrate.util import get_pointing_table
from aiapy.util.exceptions import AiapyUserWarning

__all__ = ["fix_observer_location", "update_pointing"]


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
        x=smap.meta["haex_obs"] * u.m,
        y=smap.meta["haey_obs"] * u.m,
        z=smap.meta["haez_obs"] * u.m,
        representation_type=CartesianRepresentation,
        frame=HeliocentricMeanEcliptic,
        obstime=smap.date,
    ).heliographic_stonyhurst
    # Update header
    new_meta = copy.deepcopy(smap.meta)
    new_meta["hgln_obs"] = coord.lon.to(u.degree).value
    new_meta["hglt_obs"] = coord.lat.to(u.degree).value
    new_meta["dsun_obs"] = coord.radius.to(u.m).value

    return smap._new_instance(smap.data, new_meta, plot_settings=smap.plot_settings, mask=smap.mask)


def update_pointing(smap, pointing_table=None):
    """
    Update pointing information in the `smap` header.

    This function updates the pointing information in `smap` by
    updating the ``CRPIX1, CRPIX2, CDELT1, CDELT2, CROTA2`` keywords
    in the header using the information provided in `pointing_table`.
    If `pointing_table` is not specified, the 3-hour pointing
    information is queried from the `JSOC <http://jsoc.stanford.edu/>`_.

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
    # This function can only be applied to full-resolution, full-frame images
    if not contains_full_disk(smap):
        raise ValueError("Input must be a full disk image.")
    shape_full_frame = (4096, 4096)
    if not all(d == (s * u.pixel) for d, s in zip(smap.dimensions, shape_full_frame)):
        raise ValueError(f"Input must be at the full resolution of {shape_full_frame}")
    if pointing_table is None:
        # Make range wide enough to get closest 3-hour pointing
        pointing_table = get_pointing_table(smap.date - 12 * u.h, smap.date + 12 * u.h)
    # Find row in which T_START <= T_OBS < T_STOP
    # The following notes are from a private communication with J. Serafin (LMSAL)
    # and are preserved here to explain the reasoning for selecting the particular
    # entry from the master pointing table (MPT).
    # NOTE: The 3 hour MPT entries are computed from limb fits of images with T_OBS
    # between T_START and T_START + 3hr, so any image with T_OBS equal to
    # T_START + 3hr - epsilon should still use the 3hr MPT entry for that T_START.
    # NOTE: For SDO data, T_OBS is preferred to DATE-OBS in the case of the
    # MPT, using DATE-OBS from near the slot boundary might result in selecting
    # an incorrect MPT record.
    t_obs = smap.meta.get("T_OBS")
    if t_obs is None:
        warnings.warn(
            "T_OBS key is missing from metadata. Falling back to Map.date. "
            "This may result in selecting in incorrect record from the "
            "master pointing table.",
            AiapyUserWarning,
        )
        t_obs = smap.date
    t_obs = astropy.time.Time(t_obs)
    t_obs_in_interval = np.logical_and(t_obs >= pointing_table["T_START"], t_obs < pointing_table["T_STOP"])
    if not t_obs_in_interval.any():
        raise IndexError(
            f"No valid entries for {t_obs} in pointing table "
            f'with first T_START date of {pointing_table[0]["T_START"]} '
            f'and a last T_STOP date of {pointing_table[-1]["T_STOP"]}.'
        )
    i_nearest = np.where(t_obs_in_interval)[0][0]
    w_str = f"{smap.wavelength.to(u.angstrom).value:03.0f}"
    new_meta = copy.deepcopy(smap.meta)
    # Extract new pointing parameters
    # The x0 and y0 keywords denote the location of the center
    # of the Sun in CCD pixel coordinates (0-based), but FITS WCS indexing is
    # 1-based. See Section 2.2 of
    # http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_K_AIA-SDO_FITS_Keyword_Document.pdf
    x0_mp = pointing_table[f"A_{w_str}_X0"][i_nearest].to("pix").value
    y0_mp = pointing_table[f"A_{w_str}_Y0"][i_nearest].to("pix").value
    crpix1 = x0_mp + 1
    crpix2 = y0_mp + 1
    cdelt = pointing_table[f"A_{w_str}_IMSCALE"][i_nearest].to("arcsecond / pixel").value
    # CROTA2 is the sum of INSTROT and SAT_ROT.
    # See http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_H_AIA-SDO_FITS_Keyword_Document.pdf
    # NOTE: Is the value of SAT_ROT in the header accurate?
    crota2 = pointing_table[f"A_{w_str}_INSTROT"][i_nearest] + smap.meta["SAT_ROT"] * u.degree
    crota2 = crota2.to("deg").value
    # Update headers
    for key, value in [
        ("crpix1", crpix1),
        ("crpix2", crpix2),
        ("x0_mp", x0_mp),  # x0_mp and y0_mp are not standard FITS keywords but they are
        ("y0_mp", y0_mp),  # used when respiking submaps so we update them here.
        ("cdelt1", cdelt),
        ("cdelt2", cdelt),
        ("crota2", crota2),
    ]:
        if np.isnan(value):
            # There are some entries in the pointing table returned from the JSOC that are marked as
            # MISSING. These get converted to NaNs when we cast it to an astropy quantity table. In
            # these cases, we just want to skip updating the pointing information.
            warnings.warn(
                f"Missing value in pointing table for {key}. This key will not be updated.",
                AiapyUserWarning,
            )
        else:
            new_meta[key] = value

    # sunpy map converts crota to a PCi_j matrix, so we remove it to force the
    # re-conversion.
    new_meta.pop("PC1_1")
    new_meta.pop("PC1_2")
    new_meta.pop("PC2_1")
    new_meta.pop("PC2_2")
    return smap._new_instance(smap.data, new_meta, plot_settings=smap.plot_settings, mask=smap.mask)
