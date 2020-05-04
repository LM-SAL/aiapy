import numpy as np
import astropy.units as u
from sunpy.map.sources.sdo import AIAMap
from sunpy.map.mapbase import PixelPair
import copy
from astropy.io import fits
import drms

__all__ = ['respike', 'fetch_spikes']


def respike(smap, spikes=None):
    """
    Re-insert "spikes" or "hot pixels" into level 1 AIA images

    Level 1 AIA images are, by default, "de-spiked" for hot-pixels
    to remove erroneously high intensity values, e.g. due to cosmic
    ray hits. This function re-inserts these "spikes" back into the
    image using spike location and intensity values either provided
    by the user or obtained automatically from JSOC.

    .. note:: This function should only be applied to level 1 images (i.e.
              before calling `aiapy.calibrate.register`).

    .. note:: This function modifies the `LVL_NUM`, `NSPIKES`, and `COMMENTS`
              header keywords ssuch that the resulting FITS header will differ
              from the original file.

    .. note:: If the image series of interest is large, it is advised to
              obtain the spike data via JSOC externally and specify them
              via the `spikes` keyword argument. To retrieve the coordinates
              of the positions of the spikes use the function
              `aiapy.calibrate.spike_coords`.

    Parameters
    ----------
    smap : `~sunpy.map.Map`
        A full-disk level 1 AIA image.
    spikes : array-like, with shape ``(2, N)``, optional
        Tuple of pixel positions of the spikes in the coordinate system of
        the full-disk AIA image (first row) and original intensity values. This
        can be calculated using `fetch_spikes`.

    Returns
    -------
    : `~sunpy.map.Map`
        A level 0.5 version of `mapi` with the spike data re-inserted

    See Also
    --------
    fetch_spikes
    """
    if not isinstance(smap, AIAMap):
        raise ValueError("Input must be an AIAMap.")

    if smap.meta['lvl_num'] != 1.0:
        raise ValueError('Can only apply respike procedure to level 1 data')

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
    new_data[coords.y.value.astype(int), coords.x.astype(int)] = values
    # Update metadata
    new_meta = copy.deepcopy(smap.meta)
    new_meta['lvl_num'] = 0.5
    new_meta['comments'] = f'Respike procedure applied; {values.shape[0]} hot pixels have been put back.'
    new_meta['nspikes'] = 0

    return smap._new_instance(
        new_data, new_meta, plot_settings=smap.plot_settings,)


def fetch_spikes(smap, as_coords=False):
    """
    Returns coordinates and values of removed spikes

    Returns coordinates and values of removed spikes which were removed in a
    level 1 `~sunpy.map.sources.sdo.AIAMap`. The locations of spikes are
    automatically retrieved from the JSOC.

    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
        A `~sunpy.map.Map` containing a full-disk AIA image
    as_coords : `bool`, optional
        If set True then a tuple of x and y pixel coordinates is returned.
    
    Returns
    -------
    : `~astropy.coordinates.SkyCoord` or `~sunpy.map.mapbase.PixelPair
        Locations of the removed spikes. By default, these are represented as
        pixel coordinates. If `as_coords=True`, the locations are returned in
        the projected coordinate system of the image.
    : array-like
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

    if all(d == (4096*u.pixel) for d in smap.dimensions):
        values = spikes[1, :]
        x_coords = spikes[0, :] % 4096
        y_coords = spikes[0, :] // 4096
    else:
        # NOTE: x0_mp, y0_mp are not standard FITS keywords
        xll = ((smap.meta['x0_mp'] + 1) - smap.meta['crpix1']) + 1
        yll = ((smap.meta['y0_mp'] + 1) - smap.meta['crpix2']) + 1
        _, ss_match, x_coords, y_coords = _full2part(
            spikes[0, :],
            xll,
            yll,
            smap.dimensions.x.value,
            smap.dimensions.y.value,
        )
        values = spikes[1, ss_match]

    coords = PixelPair(x_coords*u.pixel, y_coords*u.pixel)
    if as_coords:
        coords = smap.pixel_to_world(*coords)
    return coords, values


def _full2part(x, xll, yll, nx_part, ny_part):
    """
    Transform the spike locations from full-disk pixel coordinate system
    to pixel coordinate system of a cutout

    Pareameters
    -----------
    x : array-like
        1D coordinates of the spikes in full-disk pixel coordinates
    xll : array-like
        Lower-left x pixel coordinate with respect to the original full-disk map
    yll : array-like
        Lower-left y pixel coordinate
    nx_part : `int`
        length of the x dimension
    ny_part : `int`
        length of the y dimension

    Returns
    -------
    """
    x_orig = x % 4096
    y_orig = x // 4096
    x_part = x_orig - xll + 1
    y_part = y_orig - yll + 1

    ss_match = np.array(np.where(((x_part >= 0) &
                                  (y_part >= 0) &
                                  (x_part < nx_part) &
                                  (y_part < ny_part))))

    if ss_match.size == 0:
        raise IndexError('Could not find cutout pixel coordinates')

    x_part = (x_part[np.array(ss_match)]).astype(int)
    y_part = (y_part[np.array(ss_match)]).astype(int)
    out_arr = y_part*nx_part + x_part

    return out_arr, ss_match, np.squeeze(x_part), np.squeeze(y_part)
