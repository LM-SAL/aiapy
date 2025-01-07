"""
Functions for updating/fixing header keywords.
"""

import copy
import warnings

import numpy as np

import astropy.units as u

from sunpy.map import contains_full_disk

from aiapy.util.exceptions import AIApyUserWarning

__all__ = ["update_pointing"]


def update_pointing(smap, *, pointing_table):
    """
    Update the pointing information in the input map header.

    This function updates the pointing information in ``smap`` by
    updating the ``CRPIX1, CRPIX2, CDELT1, CDELT2, CROTA2`` keywords
    in the header using the information provided in ``pointing_table``.

    .. note::

        The method removes any ``PCi_j`` matrix keys in the header and
        updates the ``CROTA2`` keyword.

    .. warning::

        This function is only intended to be used for full-disk images
        at the full resolution of 4096x4096 pixels. It will raise a
        ``ValueError`` if the input map does not meet these criteria.

    Parameters
    ----------
    smap : `~sunpy.map.sources.AIAMap`
        Input map.
    pointing_table : `~astropy.table.QTable`
        Table of pointing information.
        You can get this table by calling `aiapy.calibrate.util.get_pointing_table`.

    Returns
    -------
    `~sunpy.map.sources.AIAMap`
        Updated map with pointing information.

    See Also
    --------
    `aiapy.calibrate.util.get_pointing_table`
    """
    if not contains_full_disk(smap):
        msg = "Input must be a full disk image."
        raise ValueError(msg)
    shape_full_frame = (4096, 4096)
    if not all(d == (s * u.pixel) for d, s in zip(smap.dimensions, shape_full_frame, strict=True)):
        msg = f"Input must be at the full resolution of {shape_full_frame}"
        raise ValueError(msg)
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
    # NOTE: In sunpy >=6.0, the reference_date property was introduced which, for
    # AIA maps, will always be pulled from "T_OBS"
    t_obs_in_interval = np.logical_and(
        smap.reference_date >= pointing_table["T_START"],
        smap.reference_date < pointing_table["T_STOP"],
    )
    if not t_obs_in_interval.any():
        msg = (
            f"No valid entries for {smap.reference_date} in pointing table "
            f'with first T_START date of {pointing_table[0]["T_START"]} '
            f'and a last T_STOP date of {pointing_table[-1]["T_STOP"]}.'
        )
        raise IndexError(msg)
    i_nearest = np.where(t_obs_in_interval)[0][0]
    w_str = f"{smap.wavelength.to(u.angstrom).value:03.0f}"
    new_meta = copy.deepcopy(smap.meta)
    # Extract new pointing parameters
    # The x0 and y0 keywords denote the location of the center
    # of the Sun in CCD pixel coordinates (0-based), but FITS WCS indexing is
    # 1-based. See Section 2.2 of
    # http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_K_AIA-SDO_FITS_Keyword_Document.pdf
    x0_mp = pointing_table[f"A_{w_str}_X0"][i_nearest].to_value("pix")
    y0_mp = pointing_table[f"A_{w_str}_Y0"][i_nearest].to_value("pix")
    crpix1 = x0_mp + 1
    crpix2 = y0_mp + 1
    cdelt = pointing_table[f"A_{w_str}_IMSCALE"][i_nearest].to_value("arcsecond / pixel")
    # CROTA2 is the sum of INSTROT and SAT_ROT.
    # See http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_H_AIA-SDO_FITS_Keyword_Document.pdf
    # NOTE: Is the value of SAT_ROT in the header accurate?
    crota2 = pointing_table[f"A_{w_str}_INSTROT"][i_nearest] + smap.meta["SAT_ROT"] * u.degree
    crota2 = crota2.to_value("deg")
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
                AIApyUserWarning,
                stacklevel=3,
            )
        else:
            new_meta[key] = value
    # sunpy.map.Map converts crota to a PCi_j matrix, so we remove it to force the re-conversion.
    new_meta.pop("PC1_1")
    new_meta.pop("PC1_2")
    new_meta.pop("PC2_1")
    new_meta.pop("PC2_2")
    return smap._new_instance(smap.data, new_meta, plot_settings=smap.plot_settings, mask=smap.mask)
