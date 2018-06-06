from math import cos,radians
import numpy as np
       

def compton_energy(energy,scattering_angle):
    """
    
    The scattering of photons from charged particles 
    is called Compton scattering after Arthur Compton who was the first to 
    measure photon-electron scattering in 1922. When the incoming photon 
    gives part of its energy to the electron, then the scattered photon has 
    lower energy and according to the Planck relationship has lower 
    frequency and longer wavelength. The wavelength change in such 
    scattering depends only upon the angle of scattering for a 
    given target particle
    
    Parameters:
    Inputs:
        
        energy (keV)  - input photon energy
        scattering_angle  (degrees) - compton scattering angle 
    
    Output:
        
        compton energy (keV) - location of the compton peak        
    
    
    """

#    wavelength_m = 1.0e+7 / (8.06554429e+5 * energy)
#    lamda_compton = wavelength_m + 2.43e-2 * (1. - np.cos(scattering_angle))
 #   compton_peak = 1.0e+7 / (8.06554429e+5 * lamda_compton)
    # the rest-mass energy of an electron (511 keV)
    mc2 = 511
    comp_denom = (1. + energy / mc2 *
                 (1. - np.cos(np.deg2rad(scattering_angle))))
    compton_peak = energy / comp_denom
	
    return compton_peak    





