import matplotlib.pyplot as plt
import numpy as np
from IPython.display import clear_output

def bleed(img, bit_depth):
    from numpy import zeros, floor, argsort
    shape = img.shape
    img2 = zeros(shape,dtype=float)
    img2[:,:] = img[:,:]
    
    top = 2.0**bit_depth-1
    a = argsort(img2,axis=None)
    a = a[::-1] # descending order
    j = np.array((a % shape[0]),dtype=int)
    i = np.array(floor(a/shape[0]),dtype=int)
    dval = 0.0
    
    for n in range(len(a)):
        if (img[i[n],j[n]] > top+1):
            jp1 = j[n]+1
            jm1 = j[n]-1
            val0 = img[i[n],j[n]]
            val_left = 0.5*(img[i[n],j[n]+1] + val0)
            val_right= 0.5*(img[i[n],j[n]-1] + val0)
            val1 = min((top, val_left, val_right))
            dval = val0-val1
            if (dval > 1):
                val2 = val_left/(val_left+val_right)*dval
                val3 = val_right/(val_left+val_right)*dval
                img2[i[n],j[n]] = val1
                img2[i[n],j[n]-1] += val2
                img2[i[n],j[n]+1] += val3
                return img2, dval
    return img2, dval

def bleed_wrapper(img, bit_depth, displayFreq=-1):
    result = np.array(img.shape,dtype=float)
    top = 2.0**bit_depth-1
    iternum = 0
    while (img.max() >= top):
        result, dval = bleed(img, bit_depth)
        img[:,:] = result[:,:]
        if dval ==0.0:
            return result
        if (iternum%displayFreq)==0:
            clear_output(wait=True)
            import matplotlib.pyplot as plt
            plt.imshow(np.log2(img.T), origin='lowerleft', interpolation='nearest', clim=(0,bit_depth), cmap='gnuplot2')
            plt.title('Iteration #{0:05d}: Max={1:.4g}'.format(iternum,img.max()))
            plt.show()
        iternum += 1
    return result

def read_psf(wavelnth, shape, dir='/Users/cheung/AIA/weber_psfs_rc', nocore=False):
    from astropy.io import fits 
    nx = shape[0]
    ny = shape[1]
    prefix = 'psf_'
    if nocore: 
        prefix = prefix + 'entrance_filter_only'
    #Read in PSF
    filename = '{0}/{1}_{2:04d}.fits'.format(dir,prefix,wavelnth)
    print(filename)
    hdulist = fits.open(filename)
    hdulist.info()
    psf = np.array(hdulist[0].data,dtype=float)

    #Cut out PSF and store in sparse matrix
    from scipy import sparse
    from scipy.sparse import linalg
    d_array = []
    i_array = []
    j_array = []
    count = 0
    countsub = 0 
    exclude_count = 0
    xc, yc = np.unravel_index(psf.argmax(),psf.shape)

    for ip in (np.arange(2*nx)-nx):
        clear_output(wait=True)
        for jp in (np.arange(2*ny)-ny):
            pval = psf[xc + ip, yc + jp]
            for i in np.arange(nx):
                for j in np.arange(ny):
                    if ((pval >= 3e-4) and (i+ip)>=0) and ((j+jp)>=0) and ((i+ip) < nx) and ((j+jp) < ny):
                        i_array.append(i+j*nx) 
                        j_array.append((i+ip)+(j+jp)*nx)
                        d_array.append(pval)
                        count = count+1
        print("Making sparse matrix from {0} :{1:03d} %".format(filename,int(100*(ip+nx)/(2*nx))))
    
    # C is the point spread function in matrix form
    C = sparse.csc_matrix((d_array,(i_array,j_array)),shape=(nx*ny,nx*ny),dtype='float64')
    return C, d_array, i_array, j_array

def sparse_deconv(sat, d_array, i_array, j_array, bit_depth, gap=0.7, alpha=1e-6, model=None):
    import scipy.sparse as sparse
    nx = sat.shape[0]
    ny = sat.shape[1]
    # Sparse reconstruct saturated image
    cap = float(2.0**(bit_depth-gap))
    saturated_pixels = np.where(sat > cap)
    include_pixels = np.ones(shape=[nx,ny],dtype='bool')
    include_pixels[:,:] = True
    include_pixels[saturated_pixels] = False
    exclude_pixels = np.array( (1-include_pixels) ,dtype='bool')

    list1 = ((np.arange(nx*ny,dtype='int'))[include_pixels.ravel()]).tolist()
    list2 = ((np.arange(nx*ny,dtype='int'))).tolist()

    # New matrix by excluding rows corresponding to missing pixels
    Ccsc = sparse.csc_matrix((d_array,(i_array,j_array)),shape=(nx*ny,nx*ny),dtype='float64')
    Cmod = Ccsc[list1,:][:,list2]
    print('Cmod.shape=',Cmod.shape)

    if model == None:
        from sklearn import linear_model
        #clf = linear_model.LassoCV(eps=1e-6, fit_intercept=False, positive=True)#(alpha=1e-14, max_iter=2000)
        model = linear_model.Lasso(alpha=alpha*((nx*ny)/include_pixels.sum())**2.0,
                                   max_iter=2000,fit_intercept=False, positive=True, selection='random')

    imgline = sat.ravel()
    model.fit(Cmod, imgline[include_pixels.ravel()])
    reclasso = model.coef_

    reclasso = np.reshape(reclasso, [nx,ny])
    desat = reclasso
    desat[desat < 1e-16] = 1e-16
    
    return desat, model

def pinv_deconv(sat, bit_depth, gap=0.7):
    # Sparse reconstruct saturated image
    cap = float(2.0**(bit_depth-gap))
    saturated_pixels = np.where(sat > cap)
    include_pixels = np.ones(shape=[nx,ny],dtype='bool')
    include_pixels[:,:] = True
    include_pixels[saturated_pixels] = False
    exclude_pixels = np.array( (1-include_pixels) ,dtype='bool')

    list1 = ((np.arange(nx*ny,dtype='int'))[include_pixels.ravel()]).tolist()
    list2 = ((np.arange(nx*ny,dtype='int'))).tolist()

    # New matrix by excluding rows corresponding to missing pixels
    Ccsc = sparse.csc_matrix((d_array,(i_array,j_array)),shape=(nx*ny,nx*ny),dtype='float64')
    Cmod = Ccsc[list1,:][:,list2]
    
    Cinv = np.linalg.pinv(Cmod.todense())
    deconv = Cinv.dot(sat.ravel()).reshape(nx,ny)
    
    return deconv

def compare_imgs(images,labels,types,clim=(0,17),cmap='nipy_spectral', colorbar_label='log2 [DN]',
                 ycut=[], xcut=[], bit_depth=14, savefig=None):
    numimages = len(images)
    f, axarray = plt.subplots(1,numimages, sharey=True, figsize=(2+numimages*3.5,3.5))    
    nx = images[0].shape[0]
    ny = images[0].shape[1]
    if ycut == None:
        ycut = ny/2
    mmm=[]
    for i in range(0,numimages):
        mmm.append(axarray[i].imshow(np.log2(images[i].T),origin='lowerleft',cmap=cmap,clim=clim,interpolation='nearest'))
        #ax1.plot([10,88],[52,52],color='white', linewidth=4, linestyle='-', alpha=0.54)
        #ax1.plot([10,88],[52,52],color='green', linewidth=2, linestyle='-', alpha=0.64)
        axarray[i].set_title('Sum({0})={1:4.3g}'.format(labels[i],images[i].sum()))
        axarray[i].set_xlim(0,nx)
        axarray[i].set_ylim(0,ny)
        if i == 0:
            cbar = f.colorbar(mmm[0], label=colorbar_label)

    for yy in ycut:
        axarray[0].plot([0,nx],[yy,yy],color='white', linewidth=4, linestyle='-', alpha=0.8)
        axarray[0].plot([0,nx],[yy,yy],color='black', linewidth=2, linestyle='-', alpha=0.8)
    for xx in xcut:
        axarray[0].plot([xx,xx],[0,ny],color='white', linewidth=4, linestyle='-', alpha=0.8)
        axarray[0].plot([xx,xx],[0,ny],color='black', linewidth=2, linestyle='-', alpha=0.8)
    plt.show()

    for yy in ycut:
        plt.figure(figsize=(9,3))
        for i in range(0,numimages):
            if types[i] == 'plot':
                plt.plot(images[i][:,yy],label=labels[i],linewidth=3)
            elif types[i] == 'scatter':
                plt.scatter(np.arange(nx),images[i][:,yy],label=labels[i])

        plt.plot([0,nx],[2**bit_depth,2**bit_depth],label='CCD full well',linestyle='--')
        plt.xlabel('x pixel position')
        plt.title('y = {}'.format(yy))
        plt.ylabel('DN')
        plt.yscale('log')
        plt.ylim(10,images[len(images)-1][:,yy].max()*2)
        plt.legend(loc='best')
        plt.show()

    for xx in xcut:
        plt.figure(figsize=(9,3))
        for i in range(0,numimages):
            if types[i] == 'plot':
                plt.plot(images[i][xx,:],label=labels[i],linewidth=3)
            elif types[i] == 'scatter':
                plt.scatter(np.arange(ny),images[i][xx,:],label=labels[i])

        plt.plot([0,ny],[2**bit_depth,2**bit_depth],label='CCD full well',linestyle='--')
        plt.xlabel('y pixel position')
        plt.title('x = {}'.format(xx))
        plt.ylabel('DN')
        plt.yscale('log')
        plt.ylim(10,images[len(images)-1][xx,:].max()*2)
        plt.legend(loc='best')
        plt.show()
    if len(savefig) >0:
        plt.savefig(savefig)
    plt.close()