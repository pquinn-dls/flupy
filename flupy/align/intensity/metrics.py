# -*- coding: utf-8 -*-

import numpy as np
from scipy import ndimage


EPS = np.finfo(float).eps

def sum_sq_diff(img1,img2):
    """
    
    Computes sum of the square of the differences between two imagess  
    
    Parameters
    ----------
    img1 : 2D numpy or dask array
    img2 : 2D numpy or dask array

    Returns
    -------
     sum_sq_diff : float 
    
    
    """
    
    diff = (img1-img2)
    
    return np.sum(diff**2)


def mutual_information(img1,img2,bins=None):
    """
    
    Computes mutual information between two images variate from a
    joint histogram.
    
    Parameters
    ----------
    img1 : 2D array
    img2 : 2D array
    bins : number of bins to use 
           Default = None.  If None specificed then 
           the inital estimate is set to be int(sqrt(xsize*ysize/5.))
           where xsize and ysize are the number of x any y pixels respectively        
    Returns
    -------
     mi: float  the computed similariy measure
     
     
    """

    if bins == None:
        xsize,ysize = img1.shape
        bins = int(np.sqrt(xsize*ysize/5.))

    # Convert bins counts to probability values
    #print img.shape,img2.shape,np.histogram2d(img1.ravel(),img2.ravel(),bins)
    hgram, x_edges, y_edges = np.histogram2d(img1.ravel(),img2.ravel(),bins)
    pxy = hgram / float(np.sum(hgram))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
    # Now we can do the calculation using the pxy, px_py 2D arrays
    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    return 1.0/(np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs])))

def norm_mutual_information(img1, img2, normalized=True,bins=None):
    """
    Computes (normalized) mutual information between two 1D variate from a
    joint histogram.
    Parameters
    ----------
    img1 : 2D array
        first variable
    img2 : 2D array
        second variable
    normalized - Default True
    bins : number of bins to use 
           Default = None.  If None specificed then 
           the inital estimate is set to be int(sqrt(xsize*ysize/5.))
           where xsize and ysize are the number of x any y pixels respectively
    Returns
    -------
    nmi: float
        the computed normalized similariy measure
    """
    #bins = (32, 32)
    if bins == None:
        xsize,ysize = img1.shape
        bins = int(np.sqrt(xsize*ysize/5.))

    jh = np.histogram2d(img1.ravel(), img2.ravel(), bins=bins)[0]

    # smooth the jh with a gaussian filter of given sigma
    # ndimage.gaussian_filter(jh, sigma=sigma, mode='constant',
    #                             output=jh)

    # compute marginal histograms
    jh = jh + EPS
    sh = np.sum(jh)
    jh = jh / sh
    s1 = np.sum(jh, axis=0).reshape((-1, jh.shape[0]))
    s2 = np.sum(jh, axis=1).reshape((jh.shape[1], -1))

    # Normalised Mutual Information of:
    # Studholme,  jhill & jhawkes (1998).
    # "A normalized entropy measure of 3-D medical image alignment".
    # in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.
    if normalized:
        mi = ((np.sum(s1 * np.log(s1)) + np.sum(s2 * np.log(s2)))
                / np.sum(jh * np.log(jh))) - 1
    else:
        mi = ( np.sum(jh * np.log(jh)) - np.sum(s1 * np.log(s1))
               - np.sum(s2 * np.log(s2)))

    return 1.0/mi    
    