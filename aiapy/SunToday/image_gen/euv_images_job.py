"""Program for auto-updating sdo www website images
"""

import imp
try:
    fetch = imp.load_source('fetch', 'fetch.py')
except:
    print('Failure to import from sahome.')
    fetch = imp.load_source('fetch',
                            '/Users/jinmeng/matthew/redundancy/fetch.py')
import argparse
import mov_img
from datetime import datetime, timedelta
import numpy as np
import os
from scipy import misc
from PIL import Image

waves = [94, 131, 171, 193, 211, 304, 335, 1600, 1700]
out_dir = '/opt/SunInTime/'


def make_sdo(get_fnames=False):
    """
    Function to output the most recent AIA EUV images

    Args:
        get_fnames (bool), optional: return the names of the fits files 
                                        found and not output images
    """
    logo = Image.open('/opt/SunInTime/images/large_aia_watermark_50pct_rgba.png')
    source_fnames = {}
    parser = argparse.ArgumentParser()
    parser.add_argument("--y", type=int)
    parser.add_argument("--m", type=int)
    parser.add_argument("--d", type=int)
    parser.add_argument("--h", type=int)
    parser.add_argument("--mi", type=int)
    parser.add_argument("--images", action="store_true")
    parser.add_argument("--local", action="store_true")
    args = parser.parse_args()
    now = datetime.utcnow()
    if args.y:
        now = now.replace(year=args.y)
    if args.m:
        now = now.replace(month=args.m)
    if args.d:
        now = now.replace(day=args.d)
    if args.h:
        now = now.replace(hour=args.h)
    if args.mi:
        now = now.replace(minute=args.mi)
    for wave in waves:
        onemin = timedelta(minutes=1)
        r = fetch.fetch('aia_test.lev1p5', start=now, end_or_span=30,
                        useJSOC2=True, wavelengths=wave,
                        segments='image_lev1p5', keys=['date__obs'])
        while (len(r) < 1):
            now -= onemin
            r = fetch.fetch('aia_test.lev1p5', start=now, end_or_span=30,
                            useJSOC2=True, wavelengths=wave,
                            segments='image_lev1p5', keys=['date__obs'])
        if '--local' not in sys.argv:
            r[0][-1] = 'http://jsoc.stanford.edu' + r[0][-1]
        if get_fnames:
            source_fnames[str(wave)] = r[0][-1]
        else:
            out_dir = get_out_dir(r[0][0]) # output directory
            padded_wave = str(wave).zfill(4) # zero-padded wavelength string
            date = datetime.strptime(r[0][0][:-4], '%Y-%m-%dT%H:%M:%S')
            date_str = date.strftime('%Y%m%d_%H%M%S') # date string

            im_4k = mov_img.process_img(r[0][-1]) # fullsize
            im_4k.paste(logo, (3000, 3760), mask=logo)
            im_4k.save('{}f{}.jpg'.format(out_dir, padded_wave))

            im_2k = im_4k.resize((2048, 2048), resample=Image.LANCZOS) # 2k
            #im_2k.save('{}2k_{}_{}.jpg'.format(out_dir, padded_wave, date_str))
            
            im_1k = im_4k.resize((1024, 1024), resample=Image.LANCZOS) # 1k
            im_1k.save('{}l{}.jpg'.format(out_dir, padded_wave))
            #im_1k.save('{}1k_{}_{}.jpg'.format(out_dir, padded_wave, date_str))
            img = im_1k.resize((256, 256), resample=Image.LANCZOS)
            img.save('{}t{}.jpg'.format(out_dir, padded_wave))
            if wave == 94:
                mode = 'w'
            else:
                mode = 'a'
            tf = open(out_dir + '/times.txt', mode)
            tf.write('{}: {}\n'.format(padded_wave, date_str))
            tf.close()
    r = fetch.fetch('lm_jps.m45s_nrt[$]', useJSOC2=True,
        segments='magnetogram', keys=['t_rec'])
    if '--local' not in sys.argv:
        r[0][-1] = 'http://jsoc.stanford.edu' + r[0][-1]
    img = mov_img.process_hmi(r[0][-1])
    blos = r[0][-1]
    img.save('{}f_HMImag.jpg'.format(out_dir))
    img = img.resize((1024, 1024), resample=Image.LANCZOS)
    img.save('{}l_HMImag.jpg'.format(out_dir))
    img = img.resize((256, 256), resample=Image.LANCZOS)
    img.save('{}t_HMImag.jpg'.format(out_dir))
    date = datetime.strptime(r[0][0][:-4], '%Y.%m.%d_%H:%M:%S')
    date_str = date.strftime('%Y%m%d_%H%M%S')
    tf = open(out_dir + '/times.txt', 'a')
    tf.write('blos: {}\n'.format(date_str))
    r = fetch.fetch('lm_jps.Ic_45s[$]', useJSOC2=True,
        segments='continuum', keys=['t_rec'])
    if '--local' not in sys.argv:
        r[0][-1] = 'http://jsoc.stanford.edu' + r[0][-1]
    img = mov_img.process_hmi(r[0][-1])
    img.save('{}f_HMI_cont_aiascale.jpg'.format(out_dir))
    img = img.resize((1024, 1024), resample=Image.LANCZOS)
    img.save('{}l_HMI_cont_aiascale.jpg'.format(out_dir))
    img = img.resize((256, 256), resample=Image.LANCZOS)
    img.save('{}t_HMI_cont_aiascale.jpg'.format(out_dir))
    date = datetime.strptime(r[0][0][:-4], '%Y.%m.%d_%H:%M:%S')
    date_str = date.strftime('%Y%m%d_%H%M%S')
    tf.write('cont: {}\n'.format(date_str))
    red = misc.imread(out_dir + 'f0094.jpg', mode='L')
    grn = misc.imread(out_dir + 'f0335.jpg', mode='L')
    blu = misc.imread(out_dir + 'f0193.jpg', mode='L')
    rgb = [red, grn, blu]
    img = np.stack(rgb)
    img = misc.toimage(img)
    img.paste(logo, (3000, 3760), mask=logo)
    img.save('{}f_094_335_193.jpg'.format(out_dir))
    img = img.resize((1024, 1024), resample=Image.LANCZOS)
    img.save('{}l_094_335_193.jpg'.format(out_dir))
    img = img.resize((256, 256), resample=Image.LANCZOS)
    img.save('{}t_094_335_193.jpg'.format(out_dir))
    red = misc.imread(out_dir + 'f0304.jpg', mode='L')
    grn = misc.imread(out_dir + 'f0211.jpg', mode='L')
    blu = misc.imread(out_dir + 'f0171.jpg', mode='L')
    rgb = [red, grn, blu]
    img = np.stack(rgb)
    img = misc.toimage(img)
    img.paste(logo, (3000, 3760), mask=logo)
    img.save('{}f_304_211_171.jpg'.format(out_dir))
    img = img.resize((1024, 1024), resample=Image.LANCZOS)
    img.save('{}l_304_211_171.jpg'.format(out_dir))
    img = img.resize((256, 256), resample=Image.LANCZOS)
    img.save('{}t_304_211_171.jpg'.format(out_dir))
    red = misc.imread(out_dir + 'f0211.jpg', mode='L')
    grn = misc.imread(out_dir + 'f0193.jpg', mode='L')
    blu = misc.imread(out_dir + 'f0171.jpg', mode='L')
    rgb = [red, grn, blu]
    img = np.stack(rgb)
    img = misc.toimage(img)
    img.paste(logo, (3000, 3760), mask=logo)
    img.save('{}f_211_193_171.jpg'.format(out_dir))
    img = img.resize((1024, 1024), resample=Image.LANCZOS)
    img.save('{}l_211_193_171.jpg'.format(out_dir))
    img = img.resize((256, 256), resample=Image.LANCZOS)
    img.save('{}t_211_193_171.jpg'.format(out_dir))
    blos = mov_img.process_hmi(blos, offset=0, single_channel=True)
    blos_pos = np.clip(0 - blos, 0, 127)
    blos_pos = blos_pos.astype('uint8')
    blos_neg = np.clip(blos, 0, 127)
    blos_neg = blos_neg.astype('uint8')
    gs171 = misc.imread(out_dir + 'f0171.jpg', mode='L')
    red = np.clip(0.6*blos_pos + 0.7*gs171, 0, 255).astype('uint8')
    grn = np.clip(0.2*(blos_pos + blos_neg) + 0.7*gs171,0,255).astype('uint8')
    blu = np.clip(2.0*blos_neg, 0, 255).astype('uint8')
    rgb = [red, grn, blu]
    img = np.stack(rgb)
    img = misc.toimage(img)
    img.paste(logo, (3000, 3760), mask=logo)
    img.save('{}f_HMImag_171.jpg'.format(out_dir))
    img = img.resize((1024, 1024), resample=Image.LANCZOS)
    img.save('{}l_HMImag_171.jpg'.format(out_dir))
    img = img.resize((256, 256), resample=Image.LANCZOS)
    img.save('{}t_HMImag_171.jpg'.format(out_dir))
    tf.close()
    print('Images updated. {}'.format(datetime.utcnow()))
    if get_fnames:
        return source_fnames


def get_out_dir(t):
    """
    Outputs the directory structure for the image fname
    Creates directories as needed

    Args:
        t (datetime): the date

    """
    t = out_dir + '/'.join(t[:10].split('-')) + '/'
    if not os.path.exists(t):
        os.makedirs(t)
    return t


def make_day_movs():
    """
    Makes a movie of the last 24 hours of AIA images.
    5 seconds long, 60 fps, 12.5 minutes between frames, 300 total frames

    """
    end = datetime.utcnow()
    start = end - timedelta(days=1)
    for wave in waves:
        movname = out_dir + 'daily_videos/{}.mov'.format(wave)
        r = fetch.fetch('aia_test.lev1p5', start=start,
                        end_or_span=end, useJSOC2=True,
                        wavelengths=wave, segments='image_lev1p5')
        decimate_factor = int(np.round(len(r)/300))
        r = r[::decimate_factor]
        mov_img.make_movie(r, movname=movname, downscale=(4,4))
    print('Daily movies updated. {}'.format(datetime.utcnow()))
    

if __name__ == '__main__':
    import sys
    if '--images' in sys.argv:
        make_sdo()
    elif '--movs' in sys.argv:
        make_day_movs()
