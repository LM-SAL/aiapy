"""
Create a model PSF for a particular passband.

See documentation:
    http://www.lmsal.com/sdodocs/doc?cmd=dcur&proj_num=SDOD0076&file_type=pdf
"""

import numpy as np
import astropy.units as u

__all__ = ['calc_psf']

@u.quantity_input
def getmeshinfo(wavelength: u.angstrom, use_preflightcore=None):
    """
    Returns geometric parameters for meshes in AIA filters given a passband.

    The meshpitch, or the spacing between mesh wires, is given by:

    .. math::

        d = \frac{n}{\rho} \pm \frac{n}{\rho}\frac{\delta \rho}{\rho}
        = \frac{25400 \ \mu m}{70} = 362.9 \pm 5.2 \ \mu m

    where,

    - :math: `n` : 1 inch in micrometers, which is $25400 \ \mu m$
    - :math: `\rho` : mesh density, which is 70 lines/inch
    - :math: `d` : spacing between two consecutive mesh wires

    To find the width of the wires, or meshwidth, solve for w:

    .. math::

        A = \left(1 - \frac{w}{d}\right)^2 \\
        w = 34.3 \pm 4.0 \ \mu m

    where,

    - :math: `A` : the fractional geometrical area not covered by the mesh, \
        equal to $0.82 \pm 0.02$
    - :math: `w` : the width of the mesh wires

    Parameters
    ----------
    wavelength : `str`
        Name of passband for which filter mesh info will be returned
    use_preflightcore : `int`, optional
        Type of PSF core used (the default is 0)

    Returns
    -------
    meshinfo : `dict`
        Dictionary with filter mesh info, with the following keys:

        ``"channel"``
            Channel name (`str`)
        ``"image"``
            Image used to perfrom measurements (`str`)
        ``"refimage"``
            Background image (`str`)
        ``"solarnorth"``
            Up or down, in the image (`str`)
        ``"angle1"``
            First angle (`float`)
        ``"delta1"``
            Error in first angle (`float`)
        ``"angle2"``
            Second angle (`float`)
        ``"delta2"``
            Error in second angle (`float`)
        ``"angle3"``
            Third angle (`float`)
        ``"delta3"``
            Error in third angle (`float`)
        ``"angle4"``
            Fourth angle (`float`)
        ``"delta4"``
            Error in fourth angle (`float`)
        ``"spacing"``
            Distance between diffraction spikes from entrance filter
            (`float`)
        ``"dspacing"``
            Delta spacing (`float`)
        ``"meshpitch"``
            Pitch of the mesh in micrometers (`float`)
        ``"meshwidth"``
            Width of the mesh in micrometers (`float`)
        ``"fp_spacing"``
            Distance between diffraction spikes from focal plane
            filter (`float`)
        ``"gs_width"``
            Width applied to the Gaussian such that *after* convolution
            we have the proper width (`float`)
    """

    use_preflightcore = use_preflightcore or 0  # 0 by default

    wavenames = ['94', '131', '171', '193', '211', '304', '335']
    widths_default = [4.5 for i in range(len(wavenames))]
    widths_preflight = [0.951, 1.033, 0.962, 1.512, 1.199, 1.247, 0.962]

    # Create dictionary with all values set to None.
    meshinfo = dict.fromkeys(['channel', 'image', 'refimage', 'solarnorth',
                              'angle1', 'delta1', 'angle2', 'delta2', 'angle3',
                              'delta3', 'angle4', 'delta4', 'spacing',
                              'dspacing', 'meshpitch', 'meshwidth',
                              'fp_spacing', 'gs_width'])

    # Convert to string to put into meshinfo dictionary.
    wavelength = str(wavelength)
    wavelength = wavelength.split('.')[0]

    # Fill in dictionary depending on wavelength.
    if wavelength == '211':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '211'
        meshinfo['image'] = 'AIA20101016_191038_0211.fits'
        meshinfo['refimage']='AIA20101016_190902_0211.fits'
        meshinfo['angle1'] = 49.78
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 40.08
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -40.34
        meshinfo['delta3'] = 0.02
        meshinfo['angle4'] = -49.95
        meshinfo['delta4'] = 0.02
        meshinfo['spacing'] = 19.97
        meshinfo['dspacing'] = 0.09
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.465
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 25.004795574052405
            meshinfo['fp_diffract_sum'] = 13.023499234075027

    elif wavelength == '171':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '171'
        meshinfo['image'] = 'AIA20101016_191037_0171.fits'
        meshinfo['refimage']='AIA20101016_190901_0171.fits'
        meshinfo['angle1'] = 49.81
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 39.57
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -40.13
        meshinfo['delta3'] = 0.02
        meshinfo['angle4'] = -50.38
        meshinfo['delta4'] = 0.02
        meshinfo['spacing'] = 16.26
        meshinfo['dspacing'] = 0.10
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.377
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 25.040337841745608
            meshinfo['fp_diffract_sum'] = 13.021422901836583

    elif wavelength == '193':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '193'
        meshinfo['image'] = 'AIA20101016_191056_0193.fits'
        meshinfo['refimage']='AIA20101016_190844_0193.fits'
        meshinfo['angle1'] = 49.82
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 39.57
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -40.12
        meshinfo['delta3'] = 0.03
        meshinfo['angle4'] = -50.37
        meshinfo['delta4'] = 0.04
        meshinfo['spacing'] = 18.39
        meshinfo['dspacing'] = 0.20
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.425
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 26.072038517923822
            meshinfo['fp_diffract_sum'] = 13.023190664487535

    elif wavelength == '335':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '335'
        meshinfo['image'] = 'AIA20101016_191041_0335.fits'
        meshinfo['refimage']='AIA20101016_190905_0335.fits'
        meshinfo['angle1'] = 50.40
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 39.80
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -39.64
        meshinfo['delta3'] = 0.02
        meshinfo['angle4'] = -50.25
        meshinfo['delta4'] = 0.02
        meshinfo['spacing'] = 31.83
        meshinfo['dspacing'] = 0.07
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.738
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 24.169492143169837
            meshinfo['fp_diffract_sum'] = 13.225572082618202

    elif wavelength == '304':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '304'
        meshinfo['image'] = 'AIA20101016_191021_0304.fits'
        meshinfo['refimage']='AIA20101016_190845_0304.fits'
        meshinfo['angle1'] = 49.76
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 40.18
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -40.14
        meshinfo['delta3'] = 0.02
        meshinfo['angle4'] = -49.90
        meshinfo['delta4'] = 0.02
        meshinfo['spacing'] = 28.87
        meshinfo['dspacing'] = 0.05
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.670
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 26.390290950453796
            meshinfo['fp_diffract_sum'] = 13.189189983240809

    elif wavelength == '131':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '131'
        meshinfo['image'] = 'AIA20101016_191035_0131.fits'
        meshinfo['refimage']='AIA20101016_190911_0131.fits'
        meshinfo['angle1'] = 50.27
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 40.17
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -39.70
        meshinfo['delta3'] = 0.02
        meshinfo['angle4'] = -49.95
        meshinfo['delta4'] = 0.02
        meshinfo['spacing'] = 12.37
        meshinfo['dspacing'] = 0.16
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.289
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 27.297412314098523
            meshinfo['fp_diffract_sum'] = 13.021936641040613

    elif wavelength == '94':
        ss_wave = wavenames.index(wavelength)
        meshinfo['channel'] = '94'
        meshinfo['image'] = 'AIA20101016_191039_0094.fits'
        meshinfo['refimage']='AIA20101016_190903_0094.fits'
        meshinfo['angle1'] = 49.81
        meshinfo['delta1'] = 0.02
        meshinfo['angle2'] = 40.16
        meshinfo['delta2'] = 0.02
        meshinfo['angle3'] = -40.28
        meshinfo['delta3'] = 0.02
        meshinfo['angle4'] = -49.92
        meshinfo['delta4'] = 0.02
        meshinfo['spacing'] = 8.99
        meshinfo['dspacing'] = 0.13
        meshinfo['solarnorth'] = 'UP'
        meshinfo['meshpitch'] = 363.0
        meshinfo['meshwidth'] = 34.0
        meshinfo['fp_spacing'] = 0.207
        if use_preflightcore == 1:
            meshinfo['gs_width'] = widths_preflight[ss_wave]
        else:
            meshinfo['gs_width'] = widths_default[ss_wave]
            meshinfo['en_diffract_sum'] = 25.40082267722514
            meshinfo['fp_diffract_sum'] = 13.020358384929798
    else:
        print('AIA_CALC_PSF: Input WAVELENGTH not a recognized  value. See\
        header documentation. Aborting.')
        return 0

    return meshinfo

def diffractionpattern(meshinfo, Nx, Ny, verbose=False):
    """
    Calculates PSF for diffraction and core effects given mesh parameters

    Find the intensity of a particular diffraction spike:

    .. math::

        I = \sinc^2\left(\frac{\theta \omega}{\lambda}\right) \\
        \theta = \frac{m \lambda}{d} \\
        I = \sinc^2\left(\frac{m \omega}{d}\right)

    where,

    - :math: `I` : intensity of a diffraction spike, indexed by m
    - :math: `m` : any integer (nonzero) for indexing diffraction spikes
    - :math: `w` : the width of the mesh wires
    - :math: `d` : spacing between two consecutive mesh wires
    - :math: `\lambda` : the wavelength of light

    We can model the PSF as a 2D Gaussian function of the radial distance
    r from the center, given by:

    .. math::

        I(r, \theta) = I_0 \exp\left(\frac{-r^2}{2\sigma^2}\right)

    where,

    - :math: `I_0` : the intensity of a diffraction spike, given above
    - :math: `r` : the radial distance from the center
    - :math: `\theta` : the angle around
    - :math: `\sigma` : width of Gaussians

    Parameters
    ----------
    meshinfo : `dict`
        dictionary created in getmeshinfo function that holds
        parameters for a given AIA passband

    Returns
    -------
    psfnew : `numpy.ndarray`
        2D array that is the final PSF
    """

    assert Nx % 2 == 0 # Nx should be even
    assert Ny % 2 == 0 # Ny should be even
    assert Nx >= 32 # FOV should be sufficiently large
    assert Ny >= 32 # FOV should be sufficiently large

    psf = np.zeros((Nx, Ny), dtype = float)

    # This is the width applied to the Gaussian such that *after* convolution
    # have the proper width (which is 4/3 at 1/e of max).
    width_x = meshinfo['gs_width']
    width_y = meshinfo['gs_width']
    #print('width_x = ' + str(width_x) + ' width_y = ' + str(width_y))

    # Arrays where each value is its index + 0.5.
    x = range(Nx)
    x = np.asarray(x) + 0.5
    y = range(Ny)
    y = np.asarray(y) + 0.5

    xtomult = np.ones(Nx)  # Array of 1s
    ytomult = np.ones(Ny)

    # Matrix multiplication.
    #xx = np.outer(xtomult, x)  # Transpose of line below (closer match to IDL)
    xx = np.outer(x, ytomult)
    #yy = np.outer(x, xtomult)  # Transpose of line below (closer match to IDL)
    yy = np.outer(xtomult, y)

    d = meshinfo['spacing']  # Entrance filter spacing in pixels
    dx = []  # Spacing in the x-direction, using cosine
    dy = []  # Spacing in the y-direction, using sine

    for i in range(4):

        dx.append(meshinfo['spacing']*np.cos(meshinfo[f'angle{i+1}']/180.*np.pi))
        dy.append(meshinfo['spacing']*np.sin(meshinfo[f'angle{i+1}']/180.*np.pi))

    meshratio = meshinfo['meshpitch']/meshinfo['meshwidth']
    k = 1.0/(meshratio*d)

    dpx = 0.5
    dpy = 0.5

    # This section computes the effect from entrance filter mesh diffraction.
    #j = range(-100,101,1)  # Integer for each diffraction spike
    for j in range(-100,101,1):
        tpsf = np.zeros((Nx, Ny), dtype = float)
        if j != 0:

            if j%10 ==0:
                if verbose:
                    print('First pass wide angle ' + str(j))

            intensity = (np.sinc(j*d*k))**2  # I_0 in docstring

            xc = 0.5*Nx + dx[0]*j + dpx
            yc = 0.5*Ny + dy[0]*j + dpy
            xxx = xx-xc
            yyy = yy-yc
            #psf += np.exp(-width_x*(xx-xc)**2 - width_y*(yy-yc)**2)*intensity
            tpsf += np.exp(-width_x*xxx*xxx - width_y*yyy*yyy)

            xc = 0.5*Nx + dx[1]*j + dpx
            yc = 0.5*Ny + dy[1]*j + dpy
            xxx = xx-xc
            yyy = yy-yc
            #psf += np.exp(-width_x*(xx-xc)**2 - width_y*(yy-yc)**2)*intensity
            tpsf += np.exp(-width_x*xxx*xxx - width_y*yyy*yyy)

            xc = 0.5*Nx + dx[2]*j + dpx
            yc = 0.5*Ny + dy[2]*j + dpy
            xxx = xx-xc
            yyy = yy-yc
            #psf += np.exp(-width_x*(xx-xc)**2 - width_y*(yy-yc)**2)*intensity
            tpsf += np.exp(-width_x*xxx*xxx - width_y*yyy*yyy)

            xc = 0.5*Nx + dx[3]*j + dpx
            yc = 0.5*Ny + dy[3]*j + dpy
            xxx = xx-xc
            yyy = yy-yc
            #psf += np.exp(-width_x*(xx-xc)**2 - width_y*(yy-yc)**2)*intensity
            tpsf += np.exp(-width_x*xxx*xxx - width_y*yyy*yyy)

            psf += tpsf*intensity

    if verbose:
        print('{} entrance psf.sum() {}'.format(meshinfo['channel'],psf.sum()))

    # This section applies modifications to core of the
    # entrance filter PSF and normalizes.
    #psf = psf/np.sum(psf)*0.18
    psf2 = np.exp(-width_x*(xx-0.5*Nx-dpx)**2-width_y*(yy-0.5*Ny-dpy)**2)
    #psf2 = psf2/np.sum(psf2)*0.82  # 0.82 is area not covered by mesh
    psf_entrance_filter = 0.18*psf/meshinfo['en_diffract_sum'] + 0.82*psf2/psf2.sum()

    # This section computes the effect from the
    # focal plane filter mesh diffraction.
    psf = np.zeros((Nx, Ny), dtype = float)

    d = meshinfo['fp_spacing']  # Focal plane spacing in pixels
    meshratio = meshinfo['meshpitch']/meshinfo['meshwidth']
    k = 1.0/(meshratio*d)

    dx1 = d*np.cos(45.0/180*np.pi)
    dy1 = d*np.sin(45.0/180*np.pi)

    dx2 = d*np.cos(-45.0/180*np.pi)
    dy2 = d*np.sin(-45.0/180*np.pi)

    dpx = 0.5
    dpy = 0.5

    # Loop over diffraction orders
    for j in range(1, 101, 1):
        tpsf = np.zeros((Nx, Ny), dtype = float)

        if j%10 == 0:
            if verbose:
                    print(j)

        intensity = (np.sinc(j*d*k))**2

        xc = 0.5*Nx + dx1*j + dpx
        yc = 0.5*Ny + dy1*j + dpy
        xxx = xx-xc
        yyy = yy-yc
        #psf += np.exp(-width_x*(xx-xc)**2-width_y*(yy-yc)**2)*intensity
        tpsf += np.exp(-width_x*xxx*xxx-width_y*yyy*yyy)

        xc = 0.5*Nx + dx2*j + dpx
        yc = 0.5*Ny + dy2*j + dpy
        xxx = xx-xc
        yyy = yy-yc
        #psf += np.exp(-width_x*(xx-xc)**2-width_y*(yy-yc)**2)*intensity
        tpsf += np.exp(-width_x*xxx*xxx-width_y*yyy*yyy)

        xc = 0.5*Nx - dx1*j + dpx
        yc = 0.5*Ny - dy1*j + dpy
        xxx = xx-xc
        yyy = yy-yc
        #psf += np.exp(-width_x*(xx-xc)**2-width_y*(yy-yc)**2)*intensity
        tpsf += np.exp(-width_x*xxx*xxx-width_y*yyy*yyy)

        xc = 0.5*Nx - dx2*j + dpx
        yc = 0.5*Ny - dy2*j + dpy
        xxx = xx-xc
        yyy = yy-yc
        #psf += np.exp(-width_x*(xx-xc)**2-width_y*(yy-yc)**2)*intensity
        tpsf += np.exp(-width_x*xxx*xxx-width_y*yyy*yyy)

        psf += tpsf*intensity

    if verbose:
        print('{} focal plane psf.sum() {}'.format(meshinfo['channel'],psf.sum()))

    # This section applies modifications to the core of the
    # focal plane filter PSF and normalizes.
    #psf = psf/np.sum(psf)*0.18
    psf2 = np.exp(-width_x*(xx-0.5*Nx-dpx)**2-width_y*(yy-0.5*Ny-dpy)**2)
    #psf2 = psf2/np.sum(psf2)*0.82
    psf_focal_plane = 0.18*psf/meshinfo['fp_diffract_sum'] + 0.82*psf2/psf2.sum()

    # This section generates the composite PSF and normalizes.
    psfnew = abs(np.fft.fft2(np.fft.fft2(psf_focal_plane)*np.fft.fft2(psf_entrance_filter)))

    # Shifts columns right, so last 2048 #s numbers become 1st 2048,
    # then shifts rows 2048 spaces down, so have rotated 4 quarters of array
    # halfway around.
    roll1 = np.roll(psfnew, int(0.5*Ny), axis = 1)
    roll2 = np.roll(roll1,  int(0.5*Nx), axis = 0)
    psfnew = roll2

    psfnew = psfnew/(Nx*Ny) #/np.sum(psfnew)

    return psfnew

def calc_psf(wavelength, use_preflightcore=None, Nx=4096, Ny=4096):
    """
    Calculate and return AIA PSF for a given passband as a 2D array.
    Adapted from aia_calc_psf from the SolarSoft AIA package.

    Parameters
    ----------
    wavelength : `str`
        Name of passband for which filter mesh info will be returned
    use_preflightcore : `int`, optional
        Type of PSF core used (the default is 0)
    Nx: `int`, optional
        Number of pixels in the x-dimension for output PSF image (default is 4096)
    Ny: `int`, optional
        Number of pixels in the y-dimension for output PSF image (default is 4096)

    See Also
    --------
    getmeshinfo
    diffractionpattern

    Returns
    -------
    psf : `numpy.ndarray`
        2D array that is the final PSF
    """

    use_preflightcore = use_preflightcore or 0

    meshinfo = getmeshinfo(wavelength, use_preflightcore)

    psf = diffractionpattern(meshinfo, Nx, Ny)

    return psf
