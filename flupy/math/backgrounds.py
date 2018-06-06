"""
Created on Mon Oct 08 13:31:16 2012

@author: Colin_MSI
"""

import numpy as np
import math


def poly_background_estimator(xdata, ydata, n=2, weights=None,
                         maxIterations=12, pvalue=0.9,fixed=False):
    """
    Background estimator based on orthogonal polynomials

    Input:
    xdata,ydata (numpy arrays of same length)
    pvalue  :   ratio of variance in poly to poly value at which to stop.
                    0.9 default

    Output:
            background,polynomial weights, polynomials

    S. Steenstrup J. Appl. Cryst. (1981). 14, 226--229
    "A Simple Procedure for Fitting a Background to a Certain Class of Measured Spectra"

    """
    def __generate_parameters(n, weight, xdata, ydata):
        """
        generates polynomials based on
        S. Steenstrup J. Appl. Cryst. (1981). 14, 226--229
        "A Simple Procedure for Fitting a Background to a Certain Class of Measured Spectra"
        The polynomial parameters are based on weights supplied as part of the fitting fit.
        """
    
        Npoints = xdata.size
    
        poly = np.zeros((n, Npoints), dtype=np.float64,order='F')
    
        gamma = np.zeros(n, dtype=np.float64,order='F')
        alpha = np.zeros(n, dtype=np.float64,order='F')
        beta = np.zeros(n, dtype=np.float64,order='F')
        a = np.zeros(n, dtype=np.float64,order='F')
    
        alpha[0] = (weight * xdata).sum() / weight.sum()
        beta[0] = 0.0
        poly[0] = 1.0
        poly[1] = xdata - alpha[0]
        for j in range(n):
            if j > 1:
                poly[j] = (xdata - alpha[j - 1]) * \
                    poly[j - 1] - beta[j - 1] * poly[j - 2]
    
            p = poly[j]
            g = (weight * p * p).sum()
            gamma[j] = g
            a[j] = (weight * ydata * p / g).sum()
    
            if j > 0:
                alpha[j] = (weight * xdata * p * p).sum() / g
                beta[j] = (weight * xdata * p * poly[j - 1]).sum() / gamma[j - 1]
    
        return alpha, gamma, beta, a, poly
    m = 0
    ud = False
    Npoints = xdata.size
    ydata = np.clip(ydata, 0.0001, ydata.max())
    weight = 1.0 / ydata
#     zu = np.zeros(Npoints, dtype=np.float64)
    myiter = 0
    vartest = True
    c_old = []
    if(fixed):
        maxIterations=n+1
        n=2
        myiter=2
        #
        # For n polynomials work out the coefficients...
        #
        #_alpha, gamma, _beta, c, poly = generate_parameters_hermite(n, weight, xdata, ydata)
        #_alpha, gamma, _beta, c, poly = generate_parameters(n, weight, xdata, ydata)

        #
        #
        # The sum of the polynomial components
        #
        #zu = np.arange(n)
        #index = np.arange(n)

        
    #else: 
    AIC_OLD=1.0e34            
    while myiter < maxIterations:
        myiter = myiter + 1
        #
        # For n polynomials work out the coefficients...
        #
        _alpha, gamma, _beta, c, poly = __generate_parameters(n, weight, xdata, ydata)
        #_alpha, gamma, _beta, c, poly,poly2= generate_parameters_hermite_two(n, weight, xdata, ydata)
    
        # The sum of the polynomial components
        #
        zu = (c[:, np.newaxis] * poly).sum(axis=0)
        #
        # The error estimate
        #
        y_diff = ydata - zu
        rss  = (y_diff*y_diff).sum()
        Eu = (weight * y_diff * y_diff).sum()
        #
        # Degrees of freedom
        #
        f = Npoints - n - m
        #
        # check if error has reduced enough
        #
#         AIC = 2*(Npoints + 1)+(n+1)*(math.log(2*math.pi*rss/(n+1))+1)
#         if(AIC<AIC_OLD):
#             AIC_OLD=AIC
#         else:
#             print "AIC break",n
#             break
        if (Eu < (f + m + np.sqrt(2.0 * f))):
            break
        sigma = np.sqrt(Eu / (f * gamma))
    
        m = 0.0
            #
            # update the weights...quickest way I could think of doing
            # this in numpy
            # If the y value is close to the background then the weight it is 1/background
            # if the y value is far from background then the weight is 1/(y-back)**2
            #
        index = ydata > (zu + 2.0 * np.sqrt(np.abs(zu)))
        weight = 1.0 / np.abs(zu)
        tw = ydata[index] - zu[index]
        weight[index] = 1.0 / (tw * tw)
            #weight[index] = 0.0            
        m = Npoints - index.sum()
    
            #
            # test that
            #
            # Supposed to test new and old sigma,c but they are different every time...
            #
            #c_old = c.copy()
        if(fixed):
            n=n+1
        else:    
            if myiter > 1:
                clen = min(len(c_old), len(c))
                vartest=True
                for i in range(clen):
                    if((c[i]-sigma[i]) < c_old[i] and c_old[i] < (c[i]+sigma[i])):
                        vartest=True
                    else:
                        vartest=False
                        break                                                     
    
            if vartest:
                #if abs(sigma[n - 2] / c[n - 2]) < pvalue:
#                print 'sigma',sigma, c, n, pvalue
                if abs(sigma[n-1] / c[n-1]) < pvalue:
                    if ud:
                        break
                    n = n + 1
                else:
                    if n > 3:
                        n = n - 1
                        ud = True
        c_old = c
    print "finished", index.sum(),len(zu),Npoints,n,m
    #return zu, c, poly,weight,index
    return (c[:, np.newaxis] * poly).sum(axis=0)



def LLS(data):
    """
    LLS Log sqaure root operator
    """
    return np.log(np.log(np.sqrt(data+1.0)+1.0)+1.0)

def invert_LLS(data):
    """
    Undo the LLS operator
    """
    tmp = np.exp(np.exp(data)-1.0)-1.0
    return tmp*tmp -1.0


def SNIPBack(ydata, width=15, iterations=96):
    """
    Calculates background for ydata using peak suppression.

    Args:
    -----

            ydata:      array of data of length n.
            width:      window size. (default = 5)
            iterations: Number of passes over ydata. (default = 100)

    Returns:
    --------

            Background of ydata with length n.

    """

    root2 = math.sqrt(2)
    iwidth = int(width + 0.5)

   # g_y = np.log(np.log(ydata + 1.) + 1.)
    g_y=LLS(ydata)
    g_y_copy = np.copy(g_y)
    # loop over iterations
    # for the 
    for j in range(iterations,0,-1):
        if j < 9:
            width = width/root2
            iwidth = int(width + 0.5)
            if iwidth < 1:
                iwidth = 1
        for i in range(iwidth, g_y.size - iwidth):
            mymean = (g_y[i + iwidth] + g_y[i - iwidth]) / 2.0
            if mymean < g_y[i]:
                g_y[i] = mymean
       
    high_vals = g_y_copy < g_y
    g_y[high_vals] = g_y_copy[high_vals]

    # Transform back to y data
    #return np.exp(np.exp(g_y) - 1.) - 1.
    return invert_LLS(g_y)


def remove_SNIPback(ydata,width=25, iterations=96,clip_min=0.01):
    """
    Removes background from ydata using SNIP.
    This method converts any nans to nums then removes the SNIP back
 
    Args:
    -----

            ydata:      array of data of length n.
            width:      window size. (default = 5)
            iterations: Number of passes over ydata. (default = 100)

    Returns:
    --------

            
         yback  : ydata - SNIP background 
         back   : The snip background 

    """

    ydata = np.nan_to_num(ydata)
    back  = SNIPBack(ydata, width, iterations)
    yback = ydata - back
    yback = np.clip(yback,a_min=0.01,a_max=None)
    return  yback,back

