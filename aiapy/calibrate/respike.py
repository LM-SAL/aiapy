import numpy as np
import astropy.units as u
from astropy.time import Time
from sunpy.map.sources.sdo import AIAMap
from sunpy.map import contains_full_disk


def respike(smap, spikes=None):
    """
    Processes a full-disk level 1 `~sunpy.map.sources.sdo.AIAMap` (de-spiked for hot-pixels by default) 
    into a level 1 `~sunpy.map.sources.sdo.AIAMap` with spikes re-inserted.

    .. note:: Should only be applied to level1 images.  If applied, for example to
              level 1.5 images (e.g., after calling register() to remap the image data to 0.6" per pixel
	      and rotate the FOV to true north), the result will erroneous.

	      Also this routine modifies the header information to update keywords (lvl_num=0.5, nspikes=0)
              The FITS header resulting from this procedure will therefore differ from the original
              file.


    Parameters
    ----------
    smap : `~sunpy.map.sources.sdo.AIAMap`
        A `~sunpy.map.Map` containing a full-disk AIA image

    spikes : data array (3, N_spike)
	An array with [0,:], the 1-D pixel position of and the target full-disk AIA image, 
        and [1,:] the intensity value.

    Returns
    -------
    `~sunpy.map.sources.sdo.AIAMap`:
	        A level1 copy of `~sunpy.map.sources.sdo.AIAMap` 
    and
	x_coords - x pixel coordinates of spikes 
	y_coords - y pixel coordinates of spikes

    """

    if spikes.ndim == 2:

        if not isinstance(smap, AIAMap):
            raise ValueError("Input must be an AIAMap.")

        if smap.meta['lvl_num']==1.0:
    
            if smap.meta['nspikes']==0:
                raise ValueError("No spikes were present in the level0 AIAMap...")

    ################ ADD OPTION HERE TO DOWNLOAD SPIKES FROM JSOC
            #if not isinstance(ispikedd, None):
                   # aia_jsoc_getcaldata, iindex, spikeii, spikedd, /spikes, refresh=refresh
                   #^^^ perform Fido.search to get data for AIAMap provided!!

            if ( (smap.meta['naxis1'] < 4096) or (smap.meta['naxis2'] < 4096) ):
                print("Array dimensions suggest image is a cutout. Spike coordinates will be transformed...")

                xll = ((smap.meta['x0_mp'] + 1) - smap.meta['crpix1']) + 1
                yll = ((smap.meta['y0_mp'] + 1) - smap.meta['crpix2']) + 1
                nx_part = smap.meta['naxis1']
                ny_part = smap.meta['naxis2']


                spike_coords1d_full = spikes[0,:]

                [spike_coords1d_part, ss_match, x_coords, y_coords] = full2part(spikes[0,:], 
                        xll=xll, yll=yll, nx_part=nx_part, ny_part=ny_part, mk_vec=True)

                spike_vals_full = spikes[1,:]

                if ss_match.size != 0 :
                    spike_vals_part = spike_vals_full[ss_match]
                    smap.data[y_coords, x_coords] = spike_vals_part

                return x_coords, y_coords
            else:
		#FULL DISK AIAMAP
                spike_coords1d_full = spikes[0,:]

                spike_vals_full = spikes[1,:]
                x_coords = spike_coords1d_full  %   4096
                y_coords = spike_coords1d_full  //  4096
                smap.data[y_coords, x_coords] = spike_vals_full

                smap.meta['lvl_num'] = 0.5
                smap.meta['nspikes'] = 0
                return x_coords, y_coords
        else:
            print(' Input image must be level 1 for spikes to be inserted.  Returning.')   
            return smap
    else:

        print, ' Problem with spike data array for this FITS file.  Returning.'


#------

def full2part(x, xll=None, yll=None, xsiz_orig=None, ysiz_orig=None, 
              nx_part=None, ny_part=None, mk_vec=None, **kwargs):
#function to be used to transform the spike pixel coordinates if AIAMap is a cutout

    y = kwargs.pop('y', None)
    if xsiz_orig is None : xsiz_orig = 4096
    if ysiz_orig is None : ysiz_orig = 4096
    if not y:
        x_orig = x  %   xsiz_orig
        y_orig = x  //  xsiz_orig
    else:
        x_orig = x
        y_orig = y


    x_part = x_orig - xll + 1
    y_part = y_orig - yll + 1

    ss_match = np.array(np.where( ((x_part >= 0) & (y_part >= 0) & (x_part < nx_part) & (y_part < ny_part)))) 
    n_match=len(np.squeeze(np.array(ss_match)))

    if n_match > 0 :
        x_part = (x_part[np.array(ss_match)]).astype(int)
        y_part = (y_part[np.array(ss_match)]).astype(int)
        out_arr = y_part*nx_part + x_part if mk_vec == True else np.vstack((x_part,y_part)).T
    else:
        out_arr = -1

    return out_arr, ss_match, np.squeeze(x_part), np.squeeze(y_part)


