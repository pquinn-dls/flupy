# -*- coding: utf-8 -*-
"""
Created on Thu May 03 14:49:23 2018

xrf line simulations...

@author: pq67
"""
from lmfit import Model
import numpy as np
from flupy.xray.xraylibrary2 import XrayLine,EscapeLine
import flupy.math.lineshapes as lineshapes

def generate_energy_axis_from_dict(exptdesc):
    nchannels         = exptdesc["Acquisition_instrument"]["XRF"]["Detector"]["no_of_channels"]
    gain              = exptdesc["Acquisition_instrument"]["XRF"]["Detector"]["gain"]
    offset            = exptdesc["Acquisition_instrument"]["XRF"]["Detector"]["offset"]
    energies          = generate_energy_axis(np.arange(nchannels),gain,offset)
    return energies

def getlineshapefunc(linename):
    # Start with the element function...
    func = getattr(lineshapes,linename, None)
    if func == None:
        func = getattr(lineshapes,"gaussian", None)
    return func


def generate_energy_axis(channels,gain,offset):
    if isinstance(channels, (list, tuple, np.ndarray)):
        energies          = offset + np.array(channels)*gain
    else:
        energies          = offset + np.arange(channels)*gain
    return energies


def create_xrf_lmfit_parameters(exptdesc,linedict,linefunc):
    funcmodel  = Model(linefunc)
    funcparams = funcmodel.make_params()
    # remove the area and center parameter 
    funcparams.pop("area",None)
    funcparams.pop("center",None)
    # get the function name
    linename   = linefunc.__name__
    # For each element in the linedict add a parameter
    # which will control the area
    for key in exptdesc["lineshape"]["element"][linename].keys():
        if key in funcparams:
            funcparams[key].value = \
                exptdesc["lineshape"]["element"][linename][key]["value"]
            funcparams[key].min   = \
                exptdesc["lineshape"]["element"][linename][key]["min"]
            funcparams[key].max   = \
                exptdesc["lineshape"]["element"][linename][key]["max"]
            funcparams[key].vary  = \
                exptdesc["lineshape"]["element"][linename][key]["vary"]
    for key in linedict:
        if not '+' in key:            
            funcparams.add(key,
                         **exptdesc["lineshape"]["element"][linename]["area"])
    return funcparams



def create_xrf_line(x,params,linedict,linefunc):
    funcmodel  = Model(linefunc)
    # start with the function parameters...
    funcparams = funcmodel.make_params()
    #
    # copy the function based parameters from params to funcparams
    # 
    for key in params:
        if key in funcparams:
            funcparams[key].value = params[key].value
            funcparams[key].min   = params[key].min
            funcparams[key].max   = params[key].max
            funcparams[key].vary  = params[key].vary
    result = np.zeros_like(x)  
    sigma_keys = [s for s in funcparams.keys() if 'sigma' in s]
          
    for key,linegroup in linedict.iteritems():
        if key in params:
            for line in linegroup:
                # update funcparams
                funcparams["area"].value   = params[key].value*\
                                             line.intensity
                funcparams["center"].value = line.line_energy
                for ii in sigma_keys:
                    funcparams[ii].value = np.sqrt(params[ii].value**2 + \
                             (line.line_energy*0.001))
                result  = result + funcmodel.eval(funcparams,x=x)
               
    return result


