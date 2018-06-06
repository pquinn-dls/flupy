#
#
#  General lineshapes....
#
#
import numpy as np
from scipy.special import erf,wofz,erfc
from scipy.special import gamma as gamfcn
import logging

logger = logging.getLogger(__name__)


log2 = np.log(2)
s2pi = np.sqrt(2*np.pi)
spi = np.sqrt(np.pi)
s2 = np.sqrt(2.0)
tiny = 1.e-13


def gaussian(x, area, center, sigma):
    """1 dimensional gaussian

    Parameters
    ----------
    x : array
        independent variable
    area : float
        Area of the normally distributed peak
    center : float
        center position
    sigma : float
        standard deviation
    """
    #prefactor = (area / (s2pi * sigma))
    t = (x-center)/sigma
    return np.exp(-0.5 *t*t)*(area / (s2pi * sigma))


def lorentzian(x, area, center, sigma):
    """1 dimensional lorentzian

    Parameters
    ----------
    x : array
        independent variable
    area : float
        area of lorentzian peak,
        If area is set as 1, the integral is unity.
    center : float
        center position
    sigma : float
        standard deviation
    """
    return (area / (1 + ((1.0 * x - center) / sigma) ** 2)) / (np.pi * sigma)


def lorentzian2(x, area, center, sigma):
    """1-d lorentzian squared profile

    Parameters
    ----------
    x : array
        independent variable
    area : float
        area of lorentzian squared peak,
        If area is set as 1, the integral is unity.
    center : float
        center position
    sigma : float
        standard deviation
    """

    return (area/(1 + ((x - center) / sigma)**2)**2) / (np.pi * sigma)


def voigt(x, area, center, sigma, gamma=None):
    """Convolution of gaussian and lorentzian curve.

    see http://en.wikipedia.org/wiki/Voigt_profile

    Parameters
    ----------
    x : array
        independent variable
    area : float
        area of voigt peak
    center : float
        center position
    sigma : float
        standard deviation
    gamma : float, optional
        half width at half maximum of lorentzian.
        If optional, `gamma` gets set to `sigma`
    """
    if gamma is None:
        gamma = sigma
    z = (x - center + 1j*gamma) / (sigma * s2)
    return area*wofz(z).real / (sigma*s2pi)



def pvoigt(x, area, center, sigma, fraction):
    """Linear combination  of gaussian and lorentzian

    Parameters
    ----------
    x : array
        independent variable
    area : float
        area of pvoigt peak
    center : float
        center position
    sigma : float
        standard deviation
    fraction : float
        weight for lorentzian peak in the linear combination, and (1-fraction)
        is the weight for gaussian peak.
    """
    return ((1-fraction) * gaussian(x, area, center, sigma) +
            fraction * lorentzian(x, area, center, sigma))

def asymm_pvoigt(x, area, center, sigma, fraction,skew):
    """Linear combination  of gaussian and lorentzian

    Parameters
    ----------
    x : array
        independent variable
    area : float
        area of pvoigt peak
    center : float
        center position
    sigma : float
        standard deviation
    fraction : float
        weight for lorentzian peak in the linear combination, and (1-fraction)
        is the weight for gaussian peak.
    """
    sigma_new = 2.*sigma/(1.+np.exp(skew*(x-sigma)))
    return ((1-fraction) * gaussian(x, area, center, sigma_new) +
            fraction * lorentzian(x, area, center, sigma_new))


def asymm_gaussian(x, area,center, sigma1, sigma2):
    """
    Calculate asymmetric Gaussian shape for fluorescent element line.

    Args:
    -----

            x:             array of energy values of spectrum in KeV.
            peak_position: position of characteristic line element peak
                               (xl.LineEnergy).
            sigma1:        width used for left side of peak.
            sigma2:        width used for right side of peak.

    Returns:
    --------

            Asymmetric Gaussian shape with amplitude positioned at peak_position
            with left and right hand sides of the slope at widths sigma1 and
            sigma2 respectively.

    """

    t = x - center
    t /= np.where(t < 0, sigma1, sigma2)
    g1 = np.exp(-0.5 * t * t)
    g1 /= 0.5 * (sigma1 + sigma2) * s2pi
    return g1*area
 

def split_pvoigt(x, area, center, sigma1, sigma2, fraction):
    """Split pseudo voigt - a linear combination  of gaussian and lorentzian

    Parameters
    ----------
    x : array
        independent variable
    area : float
        area of pvoigt peak
    center : float
        center position
    sigma1 : float
        standard deviation <= center position
    sigma2 : float
        standard deviation > center position
    fraction : float
        weight for lorentzian peak in the linear combination, and (1-fraction)
        is the weight for gaussian peak.
    """
    arg = (x-center)
    lor1 = (area / (1.0 + ((1.0 * arg) / sigma1) ** 2)) / (0.5*np.pi * (sigma1+sigma2))
    lor2 = (area / (1.0 + ((1.0 * arg) / sigma2) ** 2)) / (0.5*np.pi * (sigma1+sigma2))
                                                         
    prefactor = area / (s2pi * 0.5*(sigma1+sigma2))
    gauss1 = prefactor*np.exp(-0.5 *arg*arg/(sigma1*sigma1))
    gauss2 = prefactor*np.exp(-0.5 *arg*arg/(sigma2*sigma2))
    
    p1 = (1.0-fraction)*gauss1 + fraction*lor1
    p2 = (1.0-fraction)*gauss2 + fraction*lor2
    return np.where(x<=center,p1,p2)


def pearson7(x, area, center, sigma, expon):
    """Return a Pearson7 lineshape.
    Using the wikipedia definition:
    pearson7(x, center, sigma, expon) =
        amplitude*(1+arg**2)**(-expon)/(sigma*beta(expon-0.5, 0.5))
    where arg = (x-center)/sigma
    and beta() is the beta function.
    """
    arg = (x-center)/sigma
    scale = area * gamfcn(expon)/(gamfcn(0.5)*gamfcn(expon-0.5))
    return scale*(1+arg**2)**(-expon)/sigma

def EMG(x, area, center, sigma, gamma):
    """Return a exponentially modified Gaussian.
    http://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution

    area: float
       area
    center: float
       peak center
    sigma :  float
       peak standard deviation
    gamma: float
       distortion parameter
    """
    gss = gamma*sigma*sigma
#    arg1 = gamma*(center + gss/2.0 - x)
#    arg2 = (center + gss - x)/(s2*sigma)
    arg1 = gamma*(x-(center + gss/2.0))
    arg2 = (x - center - gss)/(s2*sigma)
    
    return area*(gamma/2) * np.exp(arg1) * erfc(arg2)


def skewed_voigt(x, area, center, sigma, skew,gamma=None):
    """Return a Voigt lineshape, skewed with error function.
       useful for ad-hoc Compton scatter profile
       with beta = skew/(sigma*sqrt(2))
       = voigt(x, center, sigma, gamma)*(1+erf(beta*(x-center)))

    area: float
        area of peak
    center: float
       peak center
    sigma :  float
       peak standard deviation
    skew: float
        skew < 0:  tail to low value of centroid
        skew > 0:  tail to high value of centroid
    gamma : float
        if None gamma=sigma 

    see http://en.wikipedia.org/wiki/Skew_normal_distribution

    """
    beta = skew/(s2*sigma)
    asym = 1.0 + erf(beta*(x-center))
    return asym * voigt(x, area, center, sigma,gamma)


def gausssian_step(x, area, center, sigma, peak_e):
    """
    Gauss step function is an important component in modeling compton peak.
    Use scipy erfc function. Please note erfc = 1-erf.

    Parameters
    ----------
    x : array
        data in x coordinate
    area : float
        area of gauss step function
    center : float
        center position
    sigma : float
        standard deviation
    peak_e : float
        emission energy

    Returns
    -------
    counts : array
        gaussian step peak

    References
    ----------
    .. [1]
        Rene Van Grieken, "Handbook of X-Ray Spectrometry, Second Edition,
        (Practical Spectroscopy)", CRC Press, 2 edition, pp. 182, 2007.
    """

    return (area * erfc((x - center) / (np.sqrt(2) * sigma)) /
            (2. * peak_e))


def gaussian_tail(x, area, center, sigma, gamma):
    """
    Use a gaussian tail function to simulate compton peak

    Parameters
    ----------
    x : array
        data in x coordinate
    area : float
        area of gauss tail function
        If area is set as 1, the integral is unity.
    center : float
        center position
    sigma : float
        control peak width
    gamma : float
        normalization factor

    Returns
    -------
    counts : array
        gaussian tail peak

    References
    ----------
    .. [1]
        Rene Van Grieken, "Handbook of X-Ray Spectrometry, Second Edition,
        (Practical Spectroscopy)", CRC Press, 2 edition, pp. 182, 2007.
    """

    dx_neg = np.array(x) - center
    dx_neg[dx_neg > 0] = 0

    temp_a = np.exp(dx_neg / (gamma * sigma))
    counts = (area /
              (2 * gamma * sigma * np.exp(-0.5 / (gamma**2))) * temp_a *
              erfc((x - center) / (np.sqrt(2) * sigma) +
                                 (1 / (gamma * np.sqrt(2)))))

    return counts


def compton(x, area,center,
            sigma, fwhm_corr,
            eoffset,f_step, f_tail, gamma,
            hi_f_tail, hi_gamma):
    """
    Model compton peak, which is generated as an inelastic peak and always
    stays to the left of elastic peak on the spectrum.

    Parameters
    ----------
    x : array
        energy value
    compton_amplitude : float
        area for gaussian peak, gaussian step and gaussian tail functions
    coherent_sct_energy : float
        incident energy
    fwhm_offset : float
        global fitting parameter for peak width
    compton_angle : float
        compton angle in degree
    compton_fwhm_corr : float
        correction factor on peak width
    compton_f_step : float
        weight factor of the gaussian step function
    compton_f_tail : float
        weight factor of gaussian tail on lower side
    compton_gamma : float
        normalization factor of gaussian tail on lower side
    compton_hi_f_tail : float
        weight factor of gaussian tail on higher side
    compton_hi_gamma : float
        normalization factor of gaussian tail on higher side
    epsilon : float
        energy to create a hole-electron pair
        for Ge 2.96, for Si 3.61 at 300K
        needs to double check this value

    Returns
    -------
    counts : array
        compton peak

    References
    ----------
    .. [1]
        M. Van Gysel etc, "Description of Compton peaks in energy-dispersive
        x-ray fluorescence spectra", X-Ray Spectrometry, vol. 32, pp. 139-147,
        2003.
    """

    compton_e = center+eoffset

    counts = np.zeros_like(x)

    factor = 1. / (1. + f_step + f_tail + hi_f_tail)

    value = factor * gaussian(x, area, compton_e,
                              sigma*fwhm_corr)
    counts += value

    # compton peak, step
    if f_step > 0.:
        value = factor * f_step
        value *= gausssian_step(x, area, compton_e, sigma,
                                compton_e)
        counts += value

    # compton peak, tail on the low side
    value = factor * f_tail
    value *= gaussian_tail(x, area, compton_e, sigma,
                           gamma)
    counts += value

    # compton peak, tail on the high side
    value = factor * hi_f_tail
    value *= gaussian_tail(-1. * x, area, -1. * compton_e, sigma,
                           hi_gamma)
    counts += value

    return area*counts/np.sum(counts)



def hypermet(x, area, center,
                    sigma,hyp_st_area, hyp_st_slope, 
                     hyp_lt_area, hyp_lt_slope,hyp_step_height):
    """
    hypermet function to simulate XRF peaks and/or Compton Scatter Peak
 
    Arguments
    ---------
      x          array of ordinate (energy) values
      amplitude  overall scale factor
      center     peak centroid
      sigma      Gaussian sigma
      step       step parameter for low-x erfc step [0]
      tail       amplitude of tail function         [0]
      gamma      slope of tail function             [0.1]
 
 
    Notes
    -----
    The function is given by (with error checking for_
    small values of sigma, gamma and s2 = sqrt(2) and
    s2pi = sqrt(2*pi)):
    """


    def tail(x,area,position,sigma,st_area,st_slope):
        z0 = x-position
        z1 = sigma*1.4142135623730950488
        dhelp = st_area * 0.5 * erfc((z0/z1) + 0.5 * z1/st_slope)
        dhelp = (area * dhelp)/st_slope
        f1    = 0.5 * (sigma/st_slope) * (sigma/st_slope)
        result = np.where(np.fabs(f1+z0/st_slope) <= 612, dhelp*np.exp(f1+z0/st_slope), 0.0)
        return result


    def hyp_step(x,area,position,sigma,step_height):
        z0 = x-position
        z1 = sigma*1.4142135623730950488
        result = step_height * (area/(sigma*s2pi)) * 0.5 * erfc(z0/z1)
        return result

    gauss  = gaussian(x, area, center ,sigma)
    if(hyp_st_area>0.0 and hyp_st_slope>0.0):
        short_tail = tail(x,area,center,sigma ,hyp_st_area,hyp_st_slope)
        gauss = gauss + short_tail
    if(hyp_lt_area>0.0 and hyp_lt_slope>0.0):
        long_tail = tail(x,area,center,sigma,hyp_lt_area,hyp_lt_slope)
        gauss = gauss + long_tail
    if(hyp_step_height>0.0):
        step      = hyp_step(x,area,center,sigma,hyp_step_height)
        gauss = gauss + step 
    return  gauss 



	
	