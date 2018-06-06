"""
Estimate the solid angle....
Used PyMCA method for comparison....but it doesn't account for angles etc...

"""

import numpy as np
import math


def solid_angle_fraction(detector_area=1.0,detector_distance=10.0):
    
    """
    Fraction of solid angle subtended by the detector.
    Crude calculation which doesn't account for individual sensor geometries within the detector
    so only an order of merit for multi-element detectors and a fudge factor or calibrant will be needed to
    for proper XRF quantification 	

    Parameters
    ----------
    detector_area     : float
    detector_distance : float
        
 
    Returns
    -------
    solid angle  : float
    
    """
    
    if (detector_distance > 0.0 and detector_area > 0.0):
        radius2 = detector_area/np.pi
        solidangle = 0.5 * (1.0 -  (detector_distance/math.sqrt(detector_distance**2+ radius2)))
    else:
        solidangle = 1.0
    return solidangle
