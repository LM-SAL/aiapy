"""
Image producing routine for AIA
"""

from datetime import datetime
import matplotlib.colors as colors
import numpy as np
from scipy import misc
from skimage.transform import downscale_local_mean
import os
from subprocess import Popen
from shutil import rmtree
from matplotlib.cm import get_cmap

from sunpy.map import Map
from sunpy.instr.aia import aiaprep

import sun_intensity

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw


STANDARD_INT = {'131': 6.99685, '171': 4.99803, '193': 2.9995,
                '211': 4.99801, '304': 4.99941, '335': 6.99734,
                '94': 4.99803, '1600': 2.9, '1700' : 1.0, '4500' : 0.3}
SQRT_NORM = {'131': False, '171': True, '193': False, '211': False,
             '304': False, '335': False, '94': True, '1600' : True,
            '1700' : True, '4500' : True}
MINMAX = {'131': (7, 1200),'171': (10, 6000),'193': (120, 6000),
          '211': (30, 13000), '304': (50, 2000),'335': (3.5, 1000),
          '94': (1.5, 50), '1600' : (15, 999), '1700' : (128, 5000),
          '4500' : (0, 15000)}

def process_img(fits_file, fname=None, downscale=None,
                rescale_brightness=True, side_by_side=False,
                timestamp=True, single_channel=False, 
                suppress_aia_prep=False, custom_f=lambda x: x):
    """Produces an AIA image of the Sun from a fits file

    Parameters
    ----------
    fits_file: str (Eg 'dir/image482005.fits')
        name of fits file to process
    fname : str (Eg 'dir/out_img.jpg')
        file name to save image as
        If None:
            returns a PIL image instance instead of saving directly
    downscale: tuple of two ints (Eg: (8, 8))
        downscales the data by (x, y) factor if filled
    rescale_brightness: bool
        determines if brightness correction is done
    side_by_side: bool
        make a side-by-side comparison of the scaled and not
        brightness scaled images
    timestamp: bool
        show timestamp or not
    single_channel: bool
        return image data in np array before applying the colormap
    suppress_aia_prep: bool
        not do aia prep if image below lvl 1
    custom_f: function
        custom function applied to data array
        applied early on, just after aiaprep
    """
    hdr = sun_intensity.getFitsHdr(fits_file)
    wavelength = str(hdr['wavelnth'])
    exptime = hdr['EXPTIME']
    cmap = get_cmap('sdoaia' + wavelength)
    cmap.set_bad()
    imin, imax = MINMAX[wavelength]

    themap = Map(fits_file)
    if (hdr['lvl_num'] != 1.5) and (not suppress_aia_prep):
        # perform aiaprep if data not at level 1.5
        themap = aiaprep(themap)
    data = themap.data
    data = np.flipud(data)
    data = custom_f(data) # apply custom function data
    data = data / exptime #  normalize for exposure
    norm_scale = STANDARD_INT[wavelength]
    if int(wavelength) < 999:
        dim_factor = sun_intensity.get_dim_factor(themap.date, wavelength)
    else:
        dim_factor = 1.0
    data = data * norm_scale
    if downscale:
        data = downscale_local_mean(data, downscale)

    if rescale_brightness or side_by_side:
        imin = imin / dim_factor # brightness correction
        imax = imax / dim_factor
    data[0,0] = imin # first pixel set to min
    data[0,1] = imax # second pixel sit to max

    if SQRT_NORM[wavelength]:
        norm = colors.PowerNorm(1)
        data = np.sqrt(np.clip(data, imin, imax))
    else:
        norm = colors.LogNorm(vmin=imin, vmax=imax, clip=True)
    if single_channel:
        return norm(data)
    c_img = cmap(norm(data))
    pil_img = misc.toimage(c_img[:,:,0:3]) # remove alpha channel

    width, height = pil_img.size
    if side_by_side:
        new_img = Image.new('RGB', (width * 2, height))
        new_img.paste(pil_img, (0, 0))
        second_image = process_img(fits_file, downscale=downscale,
                                   rescale_brightness=False,
                                   timestamp=False)
        new_img.paste(second_image, (width, 0))
        pil_img = new_img

    if timestamp:
        draw = ImageDraw.Draw(pil_img)
        font_height = int(height / 64)
#        font = ImageFont.truetype('/Library/Fonts/Arial.ttf',
        font = ImageFont.truetype(
              '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf',
                                  font_height)
        ts_height = { '94' : 2, '131' : 2, '171' : 2, '193' : 3, '211' : 4,
                      '304' : 3, '335' : 4, '1600' : 2, '1700' : 2 }
        tsh = ts_height[wavelength]
        draw.text((font_height, height - (tsh * font_height)),
                  'SDO/AIA- ' + wavelength + ' ' +
                  themap.date.strftime('%Y-%m-%d %H:%M:%S'),
                  font=font)

    if fname:
        pil_img.save(fname)
    else:
        return pil_img


def process_hmi(fits_file, fname=None,
                downscale=None, timestamp=True, cmap='hmimag',
                single_channel=False, offset=128, custom_f=lambda x: x):
    """Produces a HMI image of the Sun from a fits file

    Parameters
    ----------
    fits_file: str (Eg 'dir/image482005.fits')
        name of fits file to process
    rsun_obs: float
        the rsun_obs keyword
    cdelt: float
        the cdelt keyword
    fname : str (Eg 'dir/out_img.jpg')
        file name to save image as
        If None:
    downscale: tuple of two ints (Eg: (8, 8))
        downscales the data by (x, y) factor if filled
    timestamp: bool
        show timestamp or not
    cmap: str (Eg: 'Greys_r')
        colormap to use
        sunpy colormaps available to choose from
    single_channel: bool
        return image data in np array before applying the colormap
    custom_f: function
        custom function applied to data array
        applied early on, just after aiaprep
    """
    magnetgram = fits_file.find('magnetogram')
    hdr = sun_intensity.getFitsHdr(fits_file)
    themap = Map(fits_file)
    data = themap.data
    data = custom_f(data) # apply custom function data
    data = np.flipud(data)
    rsun_obs = hdr['RSUN_OBS']
    cdelt = hdr['CDELT1']
    r_pix = rsun_obs / cdelt # can subtract from this val to clip edge
    mask = sun_intensity.get_disk_mask(data.shape, r_pix)
    if magnetgram > -1:
        data = data*2 + offset
        title = 'SDO/HMI- Magnetogram'
    else:
        data = data/256
        title = 'SDO/HMI- Continuum'
    data[mask] = 0 # off disk pixel value. can be different
    if offset != 0:
        data = np.clip(data, 0, 255)
        data = data.astype('uint8')
    else:
        data = np.clip(data, -127, 127)
        data = data.astype('int8')
    if downscale:
        data = downscale_local_mean(data, downscale)

    norm = colors.SymLogNorm(1, clip=True) # what norm to try?
    cmap = get_cmap(cmap) # can try Greys or Greys_r
    if single_channel:
        return data
    pil_img = misc.toimage(data)

    width, height = pil_img.size
    if timestamp:
        draw = ImageDraw.Draw(pil_img)
        font_height = int(height / 64)
        font = ImageFont.truetype(
            '/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf',
            font_height)
        draw.text((font_height, height - (2 * font_height)),
                  title + themap.date.strftime('%Y-%m-%d %H:%M:%S'),
                  font=font, fill=255)

    if fname:
        pil_img.save(fname)
    else:
        return pil_img


def make_movie(fits_list, movname='outfile.mov', framerate=60, **kwargs):
    """Produces a movie from the list of fits files provided
    **kwargs passes args for each frame to mov_img
    """
    if not os.path.exists('/tmp/aia_movie/'):
        os.makedirs('/tmp/aia_movie/')
    for i, fits_file in enumerate(fits_list):
        fits_file = 'http://jsoc.stanford.edu' + fits_file
        process_img(fits_file,
                    fname='/tmp/aia_movie/{:05d}.jpg'.format(i),
                    **kwargs)
    pop = Popen(['ffmpeg', '-y', '-framerate', str(framerate), '-i',
             '/tmp/aia_movie/%05d.jpg', movname])
    pop.wait()
    rmtree('/tmp/aia_movie/')
