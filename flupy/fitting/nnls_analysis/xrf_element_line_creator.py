#
# imports
#
import _xraylib as xl
import numpy as np
from flupy.algorithms.xrf_calculations.lineshapes import gaussian,gaussian_tail,shelf
from flupy.algorithms.xrf_calculations.escape_line_creator import escape_peak_calculated

def characteristic_lines(x, Z, translist,transitions, sigma, tail, slope, step,\
                      detectortype,include_escape=False,correction=None):
    """
    Calculates characteristic XRF lines in given transitions by given elements Z.

    Args:
    -----

            x:              Array of energy in KeV.
            Z:              Input element number
            translist  :    index list of transitions 0=exclude, 1=include  
            transitions:    List of transitions
            sigma:          Detector precision.
            tail : Amplitude for gaussian_tail contribution.
            slope:          Slope of tail. Slope = sigma.
            step            step for hypermet
            detectortype:   Si or Ge for determining escape peak amplitude and position
            include_escape  : Include escape peak in the line profile (otherwise it can be treated seperately)
            correction      :  Attenuation correction - basically product of sample matrix, filter, detector absorption and transmission 
                            effects  

    Returns:
    --------

            Returns characteristic XRF lines for given transitions in element Z.
    """
    # the final line...
    line   = np.zeros(x.size)
    # transition lines
    tline  = np.zeros(x.size)
    # escape peak lines
    eline  = np.zeros(x.size)

    pre_attenuation_area   = np.zeros(len(translist))
    post_attenuation_area = np.zeros(len(translist))

    # Loop over the all the transitions...
    
    for i,gt in enumerate(transitions):
        tline=tline*0.0
        eline=eline*0.0
        #
        # translist says whether the transition is to be included...
        # this depends on the incident energy etc. 
        #
        if(translist[i]==0):
            pre_attenuation_area[i]=0.0
            post_attenuation_area[i]=0.0
            continue
        else:   
            for ht in gt:
                # rad rate for relative size of the transition
                amplitude = xl.RadRate(Z, ht)
                if amplitude != 0.0:
                    # energy of this transition..
                    energy = xl.LineEnergy(Z, ht)
                    # mixture of gaussian, tail + shelf
                    tline = tline + amplitude * (gaussian(x, energy, sigma)[0] +\
                           gaussian_tail(x, energy, energy, sigma, tail, slope) +\
                                         shelf(x, energy, sigma, step) )
                    # if you are including calculated escape peaks then add them in now... 
                    if(include_escape):
                        eline=eline + amplitude*escape_peak_calculated(x, energy, sigma,detectortype)[0]   
            # sum area - only for XRF peaks... 
            pre_attenuation_area[i] = pre_attenuation_area[i]+np.sum(tline)
            #
            # apply the attenutation correction after the area calculation...
            #
            if(correction!=None):
                tline=tline*correction                               
            post_attenuation_area[i] = post_attenuation_area[i]+np.sum(tline)
        # total line shape - sum of transitions lines + escape peaks..                
        line = line + tline+eline
    return line,pre_attenuation_area

