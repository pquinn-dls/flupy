# -*- coding: utf-8 -*-

import numpy as np
from lmfit import Model
from flupy.flupy.xray.compton import *

def createEnergyAxisFromDict(paramdict,lowE=0.0,highE=None):
    """
        
    Create an energy axis
    
    The dictionary input must contain:.
    paramdict["Acquisition_instrument"]["XRF"]["Detector"]["gain"]
    paramdict["Acquisition_instrument"]["XRF"]["Detector"]["offset"]
    paramdict["Acquisition_instrument"]["XRF"]["Detector"]["no_of_channels"]
    

    
    Parameters:
    Inputs:
        
        paramdict  - dictionary containing the compton function parameters
        linefuc    - the function to use for the compton peak 
    
    Output:
        
        model      - lmfit model of the compton peak with parameters initialized         
                   - lmfit model prefix set to compton_
    
    
    """
    dparams   = paramdict["Acquisition_instrument"]["XRF"]["Detector"]
    nchannels = dparams["no_of_channels"]
    gain      = dparams["gain"]
    offset    = dparams["offset"]
    
    x = np.arange(0,nchannels)
    x = x*gain + offset
    x = np.clip(x,a_min=0.0)
    return x


def updateDictWithParams(paramdict,params):
    pass


def generate_energy_axis(paramdict):
    nchannels         = paramdict["Acquisition_instrument"]["XRF"]["Detector"]["no_of_channels"]
    gain              = paramdict["Acquisition_instrument"]["XRF"]["Detector"]["gain"]
    offset            = paramdict["Acquisition_instrument"]["XRF"]["Detector"]["offset"]
    energies          = offset + np.arange(nchannels)*gain
    return energies


def trim(x, y, low, high):
    """
    Mask two arrays applying bounds to the first array.
    Parameters
    ----------
    x : array
        independent variable
    y : array
        dependent variable
    low : float
        low bound
    high : float
        high bound
    Returns
    -------
    array :
        x with new range
    array :
        y with new range
    """
    mask = (x >= low) & (x <= high)
    return x[mask], y[mask]



def define_range(data, low, high, a0, a1):
    """
    Cut x range based on offset and linear term of a linear function.
    a0 and a1 should be defined in param_dict.
    Parameters
    ----------
    data : 1D array
        raw spectrum
    low : float
        low bound in KeV
    high : float
        high bound in KeV
    a0 : float
        offset term of energy calibration
    a1 : float
        linear term of energy calibration
    Returns
    -------
    x : array
        trimmed channel number
    y : array
        trimmed spectrum according to x
    """
    x = np.arange(data.size)
    low_new = int(np.around((low - a0)/a1))
    high_new = int(np.around((high - a0)/a1))
    x0, y0 = trim(x, data, low_new, high_new)
    return x0, y0

    
def createComptonModelFromDict(paramdict,linefunc):
    """
        
    Create and initialize a comptom model
    
    lmfit parse the arguments of the linefunc to create the model parameters.

    The dictionary input must contain:.

    paramdict["Acquisition_instrument"]["beam_energy"]  (keV)
    paramdict["Acquisition_instrument"]["XRF"]["Detector"]["scattering_angle"] (degrees)
    
    A dictionary input can be used to the the parameter values.
    The method will look for the Compton parameter in paramdict["XRF"]["compton"]
    
    An example parameter is
    
    [area]
    bound_type  = Fixed, lohi or free   (dont' vary, bound, free)
    value       = 0.5  The guess value 
    min         = 0.0  Min used only if bound_type = lohi
    max         = 1.0  Max used only if bound_type = lohi

    if the arguments aren't in the dict or an empty dict is supplied
    the value default to 
    
    [param]
    bound_type  = free
    value       = 1.0  
    min         = -inf
    max         = +inf 

    
    Parameters:
    Inputs:
        
        paramdict  - dictionary containing the compton function parameters
        linefuc    - the function to use for the compton peak 
    
    Output:
        
        model      - lmfit model of the compton peak with parameters initialized         
                   - lmfit model prefix set to compton_
    
    
    """
    
    prefix="compton_"
    pdict=paramdict["XRFFit"]["compton"]
    incident_energy = paramdict["Acquisition_instrument"]["XRF"]["beam_energy"]
    angle           = paramdict["Acquisition_instrument"]["XRF"]["Detector"]["scattering_angle"]
    model = Model(linefunc,prefix=prefix)
    compton_e = compton_energy(incident_energy,angle)
    print "compton e",compton_e
    for pname in model.param_names:
        pname = model._strip_prefix(pname)
        if pname in ["beam_energy","centre","center"]:
            varyp = 'fixed'
            model.set_param_hint(pname, value=compton_e,vary=False,min=-np.inf, max=np.inf)
            continue
        if pname in pdict.keys():
            parm = pdict[pname]
            varyp= pdict[pname]["bound_type"]
            print "pname",pname,parm["value"],varyp,pdict[pname]
            if varyp =='lohi':
                model.set_param_hint(pname, value=parm["value"],vary=True,min=parm["min"], max=parm["max"])
            elif varyp=='free':
                model.set_param_hint(pname, value=parm["value"],vary=True,min=-np.inf, max=np.inf)
            else:
                model.set_param_hint(pname, value=parm["value"],vary=False,min=parm["min"], max=parm["max"])
        else:
            model.set_param_hint(pname, value=1.0,vary=True,min=-np.inf, max=np.inf)
    return model


def createScatterModelFromDict(paramdict,linefunc):
    
    """
        
    Create and initialize a scatter peak model
    
    lmfit parse the arguments of the linefunc to create the model parameters.
    
    The dictionary input must contain:.
    paramdict["Acquisition_instrument"]["beam_energy"]  (keV)
    
    The method will look for the Compton parameter in paramdict["XRF"]["Scatter"]

    An example parameter is
    
    [area]
    bound_type  = Fixed, lohi or free   (dont' vary, bound, free)
    value       = 0.5  The guess value 
    min         = 0.0  Min used only if bound_type = lohi
    max         = 1.0  Max used only if bound_type = lohi

    if the arguments aren't in the dict or an empty dict is supplied
    the value default to 
    
    [param]
    bound_type  = free
    value       = 1.0  
    min         = -inf
    max         = +inf 

    
    Parameters:
    Inputs:
        
        paramdict  - dictionary containing the compton function parameters
        linefuc    - the function to use for the scatter peak 
    
    Output:
        
        model      - lmfit model of the compton peak with parameters initialized         
                   - lmfit model prefix set to scatter_
    
    
    """

    # typical names to help identify thr peak position from the parsed linefunc
    peakcenternames = ["beam_energy","position","peak","centre","center"]
    # prefix for the lmfit model
    prefix="scatter_"
    # location of the dictionary parameter values....hopefully
    pdict=paramdict["XRFFit"]["scatter"]
    incident_energy = paramdict["Acquisition_instrument"]["XRF"]["beam_energy"]
    # create a lmfit model
    model = Model(linefunc,prefix=prefix)
    # now loop over the uninitialized model parameters and set them using the dictionary values
    for pname in model.param_names:
        pname = model._strip_prefix(pname)
        if pname in peakcenternames:
            varyp = 'fixed'
            model.set_param_hint(pname, value=incident_energy,vary=False,min=-np.inf, max=np.inf)
            continue
        if pname in pdict.keys():
            parm = pdict[pname]
            varyp= pdict[pname]["bound_type"]
            if varyp =='lohi':
                model.set_param_hint(pname, value=parm["value"],vary=True,min=parm["min"], max=parm["max"])
            elif varyp=='free':
                model.set_param_hint(pname, value=parm["value"],vary=True,min=-np.inf, max=np.inf)
            else:
                model.set_param_hint(pname, value=parm["value"],vary=False,min=parm["min"], max=parm["max"])
        else:
            # if the parameter isn't in the dictionary then just set it to 1 
            # free to move and between -inf, +inf
            model.set_param_hint(pname, value=1.0,vary=True,min=-np.inf, max=np.inf)
    return model





