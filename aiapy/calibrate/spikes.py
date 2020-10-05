import copy
import warnings

import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.wcs.utils import pixel_to_pixel
from sunpy.map.sources.sdo import AIAMap
from sunpy.map.mapbase import PixelPair
import drms

from aiapy.util import AiapyUserWarning

__all__ = ['respike', 'fetch_spikes']


def respike(smap, spikes=None):
    """
    Re-insert "spikes" or "hot pixels" into level 1 AIA images

    Level 1 AIA images are, by default, "de-spiked"
    to remove erroneously high intensity values, e.g. due to cosmic
    ray hits. This function re-inserts these "spikes" back into the
    image using spike location and intensity values either provided
    by the user or obtained automatically from
    `JSOC <http://jsoc.stanford.edu/>`_.

    .. note:: This function should only be applied to level 1 images (i.e.
              before calling `~aiapy.calibrate.register`). If the input
              image has been interpolated in any way from the original
              level 1 data, the spikes will be reinserted at the wrong
              locations.

    .. note:: This function modifies the `LVL_NUM`, `NSPIKES`, and `COMMENTS`
              header keywords ssuch that the resulting FITS header will differ
              from the original file.

    .. note:: If the image series of interest is large, it is advised to
              obtain the spike data via JSOC externally and specify them
              via the `spikes` keyword argument. To retrieve the coordinates
              of the positions of the spikes use the function
              `aiapy.calibrate.fetch_spikes`.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
        Level 1 AIA image. This can be a cutout or a full-frame image.
    spikes : array-like, with shape ``(2, N)``, optional
        Tuple of pixel positions of the spikes in the coordinate system of
        the level 1 AIA image in `smap` (first entry) and original intensity
        values (second entry). This can be calculated using `fetch_spikes`. If
        not specified, the spike positions and intensities are automatically
        queried from the JSOC.

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap`
        A level 0.5 version of `smap` with the spike data re-inserted at the
        appropriate pixels

    See Also
    --------
    fetch_spikes
    """
    if not isinstance(smap, AIAMap):
        raise ValueError("Input must be an AIAMap.")
    if smap.meta['lvl_num'] != 1.0:
        raise ValueError('Can only apply respike procedure to level 1 data')

    # Approximate check to make sure the input map has not been interpolated
    # in any way. Note that the level 1 plate scales are not exactly 0.6
    # ''/pixel, but should not differ by more than 0.1%. This is only a
    # warning because there is no exact way of determining whether an image
    # has been interpolated or not.
    nominal_scale = 0.6 * u.arcsec / u.pixel
    tol = 1e-3 * u.arcsec / u.pixel
    if not all([u.allclose(s, nominal_scale, rtol=0, atol=tol) for s in smap.scale]):
        warnings.warn(
            (f'{smap.scale} is significantly different from the expected level '
             '1 plate scale {nominal_scale}. If this map has been interpolated '
             'in any way from the level 1 image, the spike data will likely be '
             'reinserted in the incorrect pixel positions.'),
            AiapyUserWarning
        )

    # FIXME: Should raise an exception? Or just return with a warning?
    # Or better yet, why can't the logic below just handle the case of
    # no spikes?
    if smap.meta['nspikes'] == 0:
        raise ValueError('No spikes were present in the level 0 data.')

    if spikes is None:
        coords, values = fetch_spikes(smap, as_coords=False)
    else:
        coords, values = spikes

    new_data = np.copy(smap.data)
    # NOTE: the round() is needed as pixel coordinates returned by the WCS
    # transformation may be very slightly off their integer values and
    # casting them as int will sometimes result in an off-by-one error
    new_data[coords.y.value.round().astype(int),
             coords.x.value.round().astype(int)] = values
    # Update metadata
    new_meta = copy.deepcopy(smap.meta)
    new_meta['lvl_num'] = 0.5
    new_meta['comments'] = f'Respike applied; {values.shape[0]} hot pixels reinserted.'
    new_meta['nspikes'] = 0

    return smap._new_instance(
        new_data, new_meta, plot_settings=smap.plot_settings,)


def fetch_spikes(smap, as_coords=False):
    """
    Returns coordinates and values of removed spikes

    Returns coordinates and values of removed spikes which were removed in a
    level 1 AIA image. The locations of spikes are automatically retrieved
    from the JSOC.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        Level 1 AIA image. This can be a cutout or a full-disk image, but
        it should be at the original level 1 resolution.
    as_coords : `bool`, optional
        If `True`, the pixel locations are returned as a
        `~astropy.coordinates.SkyCoord` object in the projected coordinate
        system of the image.

    Returns
    -------
    `~astropy.coordinates.SkyCoord` or `~sunpy.map.mapbase.PixelPair`
        Locations of the removed spikes. By default, these are represented as
        pixel coordinates. If `as_coords=True`, the locations are returned in
        the projected coordinate system of the image.
    array-like
        Original intensity values of the spikes
    """
    if smap.wavelength not in (1600, 1700, 4500)*u.angstrom:
        series = r"aia.lev1_euv_12s"
    else:
        series = r"aia.lev1_uv_24s"
    file = drms.Client().query(
        f'{series}[{smap.date}/12s][WAVELNTH={smap.meta["wavelnth"]}]',
        seg='spikes',
    )
    _, spikes = fits.open(f'http://jsoc.stanford.edu{file["spikes"][0]}')
    spikes = spikes.data

    shape_full_frame = (4096, 4096)
    values = spikes[1, :]
    y_coords, x_coords = np.unravel_index(spikes[0, :], shape=shape_full_frame)
    # If this is a cutout, need to transform the full-frame pixel
    # coordinates into the cutout pixel coordinates and then only select
    # those in the FOV of the cutout
    if not all(d == (s*u.pixel) for d, s in zip(smap.dimensions, shape_full_frame)):
        # Construct WCS for full frame
        meta_full_frame = copy.deepcopy(smap.meta)
        meta_full_frame['crval1'] = 0.0
        meta_full_frame['crval2'] = 0.0
        # NOTE: The x0_mp and y0_mp keywords denote the location of the center
        # of the Sun in array coordinates (0-based), but FITS WCS indexing is
        # 1-based. See Section 2.2 of
        # http://jsoc.stanford.edu/~jsoc/keywords/AIA/AIA02840_K_AIA-SDO_FITS_Keyword_Document.pdf
        meta_full_frame['crpix1'] = meta_full_frame['x0_mp'] + 1
        meta_full_frame['crpix2'] = meta_full_frame['y0_mp'] + 1
        meta_full_frame['naxis1'] = shape_full_frame[0]
        meta_full_frame['naxis2'] = shape_full_frame[1]
        wcs_full_frame = smap._new_instance(
            np.zeros_like(shape_full_frame),
            meta_full_frame,
        ).wcs
        x_coords, y_coords = pixel_to_pixel(
            wcs_full_frame, smap.wcs, x_coords, y_coords)
        # Find those indices which are still in the FOV
        match = np.where(np.logical_and(
            np.logical_and(x_coords >= 0, y_coords >= 0),
            np.logical_and(x_coords < smap.dimensions.x.value,
                           y_coords < smap.dimensions.y.value)
        ))
        x_coords = x_coords[match]
        y_coords = y_coords[match]
        values = values[match]

    coords = PixelPair(x_coords*u.pixel, y_coords*u.pixel)
    if as_coords:
        coords = smap.pixel_to_world(*coords)
    return coords, values
