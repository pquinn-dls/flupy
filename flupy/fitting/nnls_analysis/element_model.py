import numpy as np
from ..general.lineshapes import gaussian,split_pvoigt
from model_tools import _gen_class_docs
from lmfit import Model
import scipy

#log2 = np.log(2)
s2pi = np.sqrt(2*np.pi)
spi = np.sqrt(np.pi)
s2 = np.sqrt(2.0)
log2 = 0.69314718055994529
sqrt2PI= np.sqrt(2.0*np.pi);
#tosigma=1.0/(2.0*sqrt(2.0*log2))

def element_peak_xrf(x, area, center,
                     delta_center, delta_sigma,
                     ratio, ratio_adjust,
                     fwhm_offset, fwhm_fanoprime,
                     e_offset, e_linear, e_quadratic,
                     epsilon=2.96):
    """
    This is a function to construct xrf element peak, which is based on
    gauss profile, but more specific requirements need to be
    considered. For instance, the standard deviation is replaced by
    global fitting parameters, and energy calibration on x is taken into
    account.

    Parameters
    ----------
    x : array
        independent variable, channel number instead of energy
    area : float
        area of gaussian function
    center : float
        center position
    delta_center : float
        adjustment to center position
    delta_sigma : float
        adjustment to standard deviation
    ratio : float
        branching ratio
    ratio_adjust : float
        value used to adjust peak height
    fwhm_offset : float
        global fitting parameter for peak width
    fwhm_fanoprime : float
        global fitting parameter for peak width
    e_offset : float
        offset of energy calibration
    e_linear : float
        linear coefficient in energy calibration
    e_quadratic : float
        quadratic coefficient in energy calibration

    Returns
    -------
    array:
        gaussian peak profile
    """
    def get_sigma(center):
        temp_val = 2 * np.sqrt(2 * np.log(2))
        return np.sqrt((fwhm_offset/temp_val)**2 + center*epsilon*fwhm_fanoprime)

    x = e_offset + x * e_linear + x**2 * e_quadratic

    return gaussian(x, area, center+delta_center,
                    delta_sigma+get_sigma(center)) * ratio * ratio_adjust


def shelf(x, area,center,sigma):
    """
    Calculates step function contribution for characteristic XRF line for
    element Z in given shell. Step function is evaluated between -inf and
    the peak position (xl.LineEnergy(Z, Shell)).

    Args:
    -----

            x:             Array of energy in KeV.
            peak_position: Energy position of peak for element Z in given shell.
            sigma:         Detector precision initially.
            S:             Size of step. Default = 0.001, initially.

    Returns:
    --------

            Returns the step function contribution for the characteristic lines
            in the fitting matrix.

    """

    #energy_shelf = c1*peak_position_c2*peak_position**2
    d_erfc = sigma * np.sqrt(2.0)
    shelfv = scipy.special.erfc((center - x) / d_erfc)
    return 0.5*shelfv*area

def truncated_exp(x, area,center,sigma, slope):
    diff = (x-center) 
    T=np.exp(diff/slope)*scipy.special.erfc((diff/(sigma*np.sqrt(2)))+ (sigma/(slope*np.sqrt(2))) )
    return T*shelf(x,1.0,center,sigma)

def hypermet_line_v2(x, area, center,
                    sigma, hyp_st_area, hyp_st_slope, hyp_step_height):
    gauss  = gaussian(x, area, center ,sigma)
    if(hyp_st_area>0.0 and hyp_st_slope>0.0):
        short_tail = truncated_exp(x,hyp_st_area,center,sigma,hyp_st_slope)
        gauss = gauss + short_tail
    if(hyp_step_height>0.0):
        step      = shelf(x,hyp_step_height,center,sigma)
        gauss     = gauss + step 
    return  gauss 
    


def tail(x,area,position,sigma,st_area,st_slope):
    z0 = x-position
    z1 = sigma*1.4142135623730950488
    dhelp = st_area * 0.5 * scipy.special.erfc((z0/z1) + 0.5 * z1/st_slope)
    dhelp = (area * dhelp)/st_slope
    f1    = 0.5 * (sigma/st_slope) * (sigma/st_slope)
    result = np.where(np.fabs(f1+z0/st_slope) <= 612, dhelp*np.exp(f1+z0/st_slope), 0.0)
    return result


def hyp_step(x,area,position,sigma,step_height):

    z0 = x-position
    z1 = sigma*1.4142135623730950488
    result = step_height * (area/(sigma*sqrt2PI)) * 0.5 * scipy.special.erfc(z0/z1)
    return result

def hypermet_int(x, area, center,delta_center, delta_sigma,
                     ratio, ratio_adjust,
                     fwhm_offset, fwhm_fanoprime,
                     e_offset, e_linear, e_quadratic, hyp_st_area, hyp_st_slope, hyp_lt_area, hyp_lt_slope,hyp_step_height,epsilon=2.96):
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
    def get_sigma(center):
        temp_val = 2 * np.sqrt(2 * np.log(2))
        return np.sqrt((fwhm_offset/temp_val)**2 + center*epsilon*fwhm_fanoprime)

    sigma = get_sigma(center)
    sigma = sigma + delta_sigma
    center = center + delta_center
 
    z0 = x-center
    z1 = sigma*1.4142135623730950488
	
    x = e_offset + x * e_linear + x**2 * e_quadratic
    result = np.zeros_like(x)

    #gaus  = gaussian(x, area, center+delta_center,delta_sigma+sigma)

    #short_tail = tail(x,area,center + delta_center,sigma + delta_sigma ,hyp_st_area,hyp_st_slope)
    
    #long_tail = tail(x,area,center + delta_center,sigma + delta_sigma ,hyp_lt_area,hyp_lt_slope)

    gaus = ((area / (s2pi * sigma)) * np.exp(-1* z0 ** 2 / (2. * sigma ** 2)))

    # long tail
    long_tail = np.zeros_like(x)
    if ((hyp_lt_slope != 0) and (hyp_lt_area != 0)):
        dhelp = hyp_lt_area * 0.5 * scipy.special.erfc((z0/z1) + 0.5 * z1/hyp_lt_slope)
        dhelp = (area * dhelp)/hyp_lt_slope
        f1    = 0.5 * (sigma/hyp_lt_slope) * (sigma/hyp_lt_slope)
        long_tail = np.where(np.fabs(f1+z0/hyp_lt_slope) <= 612, dhelp*np.exp(f1+z0/hyp_lt_slope), 0.0)
	gaus = gaus + long_tail

    # short tail
    short_tail = np.zeros_like(x)

    if ((hyp_st_slope != 0) and (hyp_st_area != 0)):
        dhelp = hyp_st_area * 0.5 * scipy.special.erfc((z0/z1) + 0.5 * z1/hyp_st_slope)
        dhelp = (area * dhelp)/hyp_st_slope
        f1    = 0.5 * (sigma/hyp_st_slope) * (sigma/hyp_st_slope)
        short_tail = np.where(np.fabs(f1+z0/hyp_st_slope) <= 612, dhelp*np.exp(f1+z0/hyp_st_slope), 0.0)
	gaus = gaus + short_tail


    # step
    step = np.zeros_like(x)

    if (hyp_step_height !=0 and sigma !=0.0):
        step= hyp_step_height * (area/(sigma*sqrt2PI)) * 0.5 * scipy.special.erfc(z0/z1)
	gaus=gaus + step

    return gaus * ratio * ratio_adjust


def hypermet(x, area, center,delta_center, delta_sigma,
                     ratio, ratio_adjust,
                     fwhm_offset, fwhm_fanoprime,
                     e_offset, e_linear, e_quadratic, hyp_st_area, hyp_st_slope, hyp_lt_area, hyp_lt_slope,hyp_step_height,epsilon=2.96):
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
    def get_sigma(center):
        temp_val = 2 * np.sqrt(2 * np.log(2))
        return np.sqrt((fwhm_offset/temp_val)**2 + center*epsilon*fwhm_fanoprime)
 
    sigma = get_sigma(center)

    x = e_offset + x * e_linear + x**2 * e_quadratic

    gaus  = gaussian(x, area, center+delta_center,delta_sigma+sigma)

    short_tail = tail(x,area,center + delta_center,sigma + delta_sigma ,hyp_st_area,hyp_st_slope)
    
    long_tail = tail(x,area,center + delta_center,sigma + delta_sigma ,hyp_lt_area,hyp_lt_slope)

    step      = hyp_step(x,area,center+delta_center,sigma+delta_sigma,hyp_step_height)

    return (gaus + step + short_tail+long_tail) * ratio * ratio_adjust


def xrf_line(x, area,
                     fwhm_offset, fwhm_fanoprime,
                     e_offset, e_linear, e_quadratic, 
                      hyp_st_area, hyp_st_slope, 
                     hyp_lt_area, hyp_lt_slope,hyp_step_height,
                     epsilon=2.96,emission_lines=None):
    #
    # lmfit parameters and models spend quite a bit of time
    # stripping strigns and intrepting  expr constraints.
    # To reduce the overhead we produce a single xrf line per element
    # rather than create a spectra from 5-10 models of individual lines. 
    # emission lines = list of (line,reverse_line_dict[line],lineE,cs))  which is all you need to work out 
    # the relative peak positions and scalings 
    #
    # The core parameters for the line are the hypermet parameters 
    # The keywords pass the element atomic number, lines and cross_sections of the lines.
    # 
    #
    temp_val = 2 * np.sqrt(2 * np.log(2))
    sig_loc = (fwhm_offset/temp_val)**2
    # calculate the x values
    x = e_offset + x * e_linear + x**2 * e_quadratic
    # create an empty line
    total_line = np.zeros_like(x)
    
    for index,item in enumerate(emission_lines):

        # reference area with respect to the first line which should be the main line
        # e.g. for Fe K this would be Ka1
        ratio_v = item[3]/emission_lines[0][3]
        # energy of line
        energy  = item[2]
        sigma = np.sqrt(sig_loc + energy*epsilon*fwhm_fanoprime)
        total_line=total_line +hypermet_line(x, area*ratio_v, energy,
                    sigma,hyp_st_area, hyp_st_slope, 
                     hyp_lt_area, hyp_lt_slope,hyp_step_height)

    return total_line
        

def xrf_line_v2(x, area,
                     fwhm_offset, fwhm_fanoprime,
                     e_offset, e_linear, e_quadratic, 
                      hyp_st_area, hyp_st_slope, 
                     hyp_step_height,
                     epsilon=2.96,emission_lines=None):
    #
    # lmfit parameters and models spend quite a bit of time
    # stripping strigns and intrepting  expr constraints.
    # To reduce the overhead we produce a single xrf line per element
    # rather than create a spectra from 5-10 models of individual lines. 
    # emission lines = list of (line,reverse_line_dict[line],lineE,cs))  which is all you need to work out 
    # the relative peak positions and scalings 
    #
    # The core parameters for the line are the hypermet parameters 
    # The keywords pass the element atomic number, lines and cross_sections of the lines.
    # 
    #
    temp_val = 2 * np.sqrt(2 * np.log(2))
    sig_loc = (fwhm_offset/temp_val)**2
    # calculate the x values
    x = e_offset + x * e_linear + x**2 * e_quadratic
    # create an empty line
    total_line = np.zeros_like(x)
    
    for index,item in enumerate(emission_lines):

        # reference area with respect to the first line which should be the main line
        # e.g. for Fe K this would be Ka1
        ratio_v = item[3]/emission_lines[0][3]
        # energy of line
        energy  = item[2]
        sigma = np.sqrt(sig_loc + energy*epsilon*fwhm_fanoprime)
        total_line=total_line +hypermet_line_v2(x, area*ratio_v, energy,
                    sigma, hyp_st_area, hyp_st_slope, hyp_step_height)

    return total_line
        
def xrf_line_v3(x, area,
                      fwhm_fanoprime,
                     e_offset, e_linear, e_quadratic, 
                      sigma1,sigma2,fraction,
                     epsilon=2.96,emission_lines=None):
    #
    # lmfit parameters and models spend quite a bit of time
    # stripping strigns and intrepting  expr constraints.
    # To reduce the overhead we produce a single xrf line per element
    # rather than create a spectra from 5-10 models of individual lines. 
    # emission lines = list of (line,reverse_line_dict[line],lineE,cs))  which is all you need to work out 
    # the relative peak positions and scalings 
    #
    # The core parameters for the line are the hypermet parameters 
    # The keywords pass the element atomic number, lines and cross_sections of the lines.
    # 
    #
    temp_val = 2 * np.sqrt(2 * np.log(2))
    #sig_loc = (fwhm_offset/temp_val)**2
    # calculate the x values
    x = e_offset + x * e_linear + x**2 * e_quadratic
    # create an empty line
    total_line = np.zeros_like(x)
    #sigma1 = 
    for index,item in enumerate(emission_lines):

        # reference area with respect to the first line which should be the main line
        # e.g. for Fe K this would be Ka1
        ratio_v = item[3]/emission_lines[0][3]
        # energy of line
        energy  = item[2]
        sigma1a = np.sqrt(sigma1**2 + energy*epsilon*fwhm_fanoprime)
        sigma2a = np.sqrt(sigma2**2 + energy*epsilon*fwhm_fanoprime)

        total_line=total_line +split_pvoigt(x, area*ratio_v, energy,
                    sigma1a,sigma2a,fraction)

    return total_line
        

    
    
def hypermet_line(x, area, center,
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



class ElementModel(Model):

    __doc__ = _gen_class_docs(element_peak_xrf)

    def __init__(self, *args, **kwargs):
#        super(ElementModel, self).__init__(element_peak_xrf, *args, **kwargs)
        super(ElementModel, self).__init__(xrf_line_v3, *args, **kwargs)

        self.set_param_hint('epsilon', value=2.96, vary=False)

