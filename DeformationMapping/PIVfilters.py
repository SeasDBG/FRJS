#Python implimentation of PIV filters from openPIV
#import replace_nan
#from replace_nan import replace_nans
import numpy as np
from scipy.signal import convolve
import cython  
    
def _gaussian_kernel( half_width=1 ):
    """A normalized 2D Gaussian kernel array
    
    Parameters
    ----------
    half_width : int
        the half width of the kernel. Kernel
        has shape 2*half_width + 1 (default half_width = 1, i.e. 
        a Gaussian of 3 x 3 kernel)
        
    Examples
    --------
    
    >>> from openpiv.filters import _gaussian_kernel
    >>> _gaussian_kernel(1)
    array([[ 0.04491922,  0.12210311,  0.04491922],
       [ 0.12210311,  0.33191066,  0.12210311],
       [ 0.04491922,  0.12210311,  0.04491922]])
   
    
    """
    size = int(half_width)
    x, y = np.mgrid[-half_width:half_width+1, -half_width:half_width+1]
    g = np.exp(-(x**2/float(half_width)+y**2/float(half_width)))
    return g / g.sum()

def gaussian_kernel(sigma, truncate=4.0):
    """
    Return Gaussian that truncates at the given number of standard deviations. 
    """

    sigma = float(sigma)
    radius = int(truncate * sigma + 0.5)

    x, y = np.mgrid[-radius:radius+1, -radius:radius+1]
    sigma = sigma**2

    k = 2*np.exp(-0.5 * (x**2 + y**2) / sigma)
    k = k / np.sum(k)

    return k


def gaussian( u, v, half_width=1) :
    """Smooths the velocity field with a Gaussian kernel.
    
    Parameters
    ----------
    u : 2d np.ndarray
        the u velocity component field
        
    v : 2d np.ndarray
        the v velocity component field
        
    half_width : int
        the half width of the kernel. Kernel
        has shape 2*half_width+1, default = 1
        
    Returns
    -------
    uf : 2d np.ndarray
        the smoothed u velocity component field
        
    vf : 2d np.ndarray
        the smoothed v velocity component field    
        
    """
    g = _gaussian_kernel( half_width=half_width )
    uf = convolve( u, g, mode='same')
    vf = convolve( v, g, mode='same')
    return uf, vf
    
DTYPEf = np.float
#typedef np.float_t DTYPEf_t
DTYPEi = np.int
#typedef np.int_t DTYPEi_t

def replace_nans(array, max_iter, tol, kernel_size=2, method='disk'):
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=DTYPEf ) 
    
    inans = np.empty([array.shape[0]*array.shape[1]], dtype=DTYPEi)
    jnans = np.empty([array.shape[0]*array.shape[1]], dtype=DTYPEi)

    iter_seeds = np.zeros(max_iter, dtype=DTYPEi)

    # indices where array is NaN
    inans, jnans = [x.astype(DTYPEi) for x in np.nonzero(np.isnan(array))]

    # number of NaN elements
    n_nans = len(inans)
    
    # arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=DTYPEf)
    replaced_old = np.zeros( n_nans, dtype=DTYPEf)
    
    # depending on kernel type, fill kernel array
    if method == 'localmean':
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1.0

    elif method == 'disk':
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                if ((kernel_size-i)**2 + (kernel_size-j)**2)**0.5 <= kernel_size:
                    kernel[i,j] = 1.0
                else:
                    kernel[i,j] = 0.0

    elif method == 'distance': 
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                if ((kernel_size-i)**2 + (kernel_size-j)**2)**0.5 <= kernel_size:
                    kernel[i,j] = -1*(((kernel_size-i)**2 + (kernel_size-j)**2)**0.5 - ((kernel_size)**2  (kernel_size)**2)**0.5)
                else:
                    kernel[i,j] = 0.0
    else:
        raise ValueError( 'method not valid. Should be one of `localmean`, `disk` or `distance`.')
        
    # make several passes
    # until we reach convergence 
    for it in range(max_iter):
    
        # for each NaN element
        for k in range(n_nans):
            i = inans[k]
            j = jnans[k]
           
            #init to 0.0
            replaced_new[k] = 0.0
            n = 0.0
            
            # loop over the kernel
            for I in range(2*kernel_size+1):
                for J in range(2*kernel_size+1):
                   
                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:
                                                
                            # if the neighbour element is not NaN itself.
                            if not np.isnan(array[i+I-kernel_size, j+J-kernel_size]):

                                # do not bother with 0 kernel values
                                if kernel[I, J] != 0:

                                    # convolve kernel with original array
                                    replaced_new[k] = replaced_new[k] + array[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + kernel[I,J]

            
            # divide value by effective number of added elements
            if n > 0:
                replaced_new[k] = replaced_new[k] / n
            else:
                replaced_new[k] = np.nan

        # bulk replace all new values in array
        for k in range(n_nans):
            array[inans[k],jnans[k]] = replaced_new[k]

        # check if mean square difference between values of replaced 
        #elements is below a certain tolerance
        if np.mean((replaced_new-replaced_old)**2) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return array
    
def replace_outliers( u, v, method='localmean', max_iter=5, tol=1e-3, kernel_size=1):
    uf = u.copy()
    vf = v.copy()
    uf = replace_nans( uf, method=method, max_iter=max_iter, tol=tol, kernel_size=kernel_size )
    vf = replace_nans( vf, method=method, max_iter=max_iter, tol=tol, kernel_size=kernel_size )
    
    return uf, vf
