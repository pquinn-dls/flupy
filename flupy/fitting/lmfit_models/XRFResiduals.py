# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model,Parameters,minimize
import numpy as np
from flupy.xray.xraylibrary2 import XrayLine,EscapeLine,ComptonLine,ElasticLine
from copy import copy
import inspect


def simple_parse_params(func,independent_vars,prefix=None):
    """Build parameters from function arguments."""
    if func is None:
        return
    if hasattr(func, 'argnames') and hasattr(func, 'kwargs'):
        pos_args = func.argnames[:]
        kw_args = {}
        for name, defval in func.kwargs:
            kw_args[name] = defval
    else:
        try:  # PY3
            argspec = inspect.getfullargspec(func)
        except AttributeError:  # PY2
            argspec = inspect.getargspec(func)

        pos_args = argspec.args
        kw_args = {}
        if argspec.defaults is not None:
            for val in reversed(argspec.defaults):
                kw_args[pos_args.pop()] = val

    # default independent_var = 1st argument
    if independent_vars is None:
        independent_vars = [pos_args[0]]

    # default param names: all positional args
    # except independent variables
    _param_root_names=None
    if _param_root_names is None:
        _param_root_names = pos_args[:]
        for p in independent_vars:
            if p in _param_root_names:
                _param_root_names.remove(p)

    used_names = []
    names = []    
    if prefix is None:
        prefix = ''
    for pname in _param_root_names:
        used_names.append("%s%s" % (prefix, pname))
        names.append("%s" % (pname))
    
    return names,used_names,func.__name__
    


def build_xrf_params(exptdesc,linedict,xrffunc,comptonfunc,elasticfunc):
    params = Parameters()

    xrfnames,xrfnames_used,xrffuncname =\
    simple_parse_params(xrffunc,["x","center","area"],prefix="xrf_")
    avail_hints = exptdesc["lineshape"]["element"][xrffuncname].keys()
    for name,pname in zip(xrfnames,xrfnames_used):
        if name in avail_hints:
            params.add(pname, **exptdesc["lineshape"]\
                                      ["element"][xrffuncname][name])
        else:
            params.add(pname,value = 1.0)

    comptonnames,comptonnames_used,comptonfuncname=\
     simple_parse_params(comptonfunc,["x","center","area"],prefix="compton_")
    avail_hints = exptdesc["lineshape"]["compton"][comptonfuncname].keys()
    for name,pname in zip(comptonnames,comptonnames_used):
        if name in avail_hints:
            params.add(pname, **exptdesc["lineshape"]\
                                      ["compton"][comptonfuncname][name])
        else:
            params.add(pname,value = 1.0)

    elasticnames,elasticnames_used,elasticfuncname =\
    simple_parse_params(elasticfunc,["x","center","area"],prefix="elastic_")
    avail_hints = exptdesc["lineshape"]["elastic"][elasticfuncname].keys()
    for name,pname in zip(elasticnames,elasticnames_used):
        if name in avail_hints:
            params.add(pname, **exptdesc["lineshape"]\
                                      ["elastic"][elasticfuncname][name])
        else:
            params.add(pname,value = 1.0)
   

    # exclude pileups from the parameters
    for key in linedict:
        if not '+' in key:  
            if "Elastic" in key:
                params.add(key,
                    **exptdesc["lineshape"]["elastic"]\
                    [elasticfuncname]["area"])
            elif "Compton" in key:
                params.add(key,
                    **exptdesc["lineshape"]["compton"]\
                    [comptonfuncname]["area"])
            else:
                params.add(key,
                    **exptdesc["lineshape"]["element"]\
                    [xrffuncname]["area"])

  #  params.add("pileup_factor",value=0.05,min=1.0e-5,max=2.,vary=True)
 #   params.add("fano",**exptdesc["xrf_fitting"]["params"]["fano"])

    customarg={}
    customarg["xrf"]={}
    customarg["xrf"]["function"]=xrffunc
    customarg["xrf"]["funcparams_orig"]=xrfnames
    customarg["xrf"]["funcparams_used"]=xrfnames_used
    
    customarg["elastic"]={}
    customarg["elastic"]["function"]=elasticfunc
    customarg["elastic"]["funcparams_orig"]  = elasticnames
    customarg["elastic"]["funcparams_used"] = elasticnames_used
    
    customarg["compton"]={}
    customarg["compton"]["function"]=comptonfunc
    customarg["compton"]["funcparams_orig"]=comptonnames
    customarg["compton"]["funcparams_used"]=comptonnames_used

    return params,customarg

#def generate_params

def xrfcompositefunc(params,x,data,customarg,xrffunc,comptonfunc,elasticfunc,linedict):
    #
    # Developed to minimize
    # the parameter handling overhead of using Models for 
    # multi-model fitting
    #
    # xrf_fitted_param_names
    # elastic_fitted_param_names
    # compton_fitted_param_names
    #
    #xrffunc     = customarg["xrf"]["function"]
    xrffunc_po  = customarg["xrf"]["funcparams_orig"]
    xrffunc_pu  = customarg["xrf"]["funcparams_used"]
    xrfparams = {o:params[u].value  for u,o in zip(xrffunc_pu,xrffunc_po)}
    

    #elasticfunc = customarg["elastic"]["function"]
    elastic_po  = customarg["elastic"]["funcparams_orig"]
    elastic_pu  = customarg["elastic"]["funcparams_used"]
    elasticparams = {o:params[u].value  for u,o in zip(elastic_pu,elastic_po)}

    #comptonfunc = customarg["compton"]["function"]
    compton_po  = customarg["compton"]["funcparams_orig"]
    compton_pu  = customarg["compton"]["funcparams_used"]
    comptonparams = {o:params[u].value  for u,o in zip(compton_pu,compton_po)}

    
#    xrf_sigma_keys     = [s for s in xrfparams if 'sigma' in s]
#    elastic_sigma_keys = [s for s in elasticparams if 'sigma' in s]
#    compton_sigma_keys = [s for s in comptonparams if 'sigma' in s]
    
#    pileup_factor = params["pileup_factor"].value 
 #   fano          = params["fano"].value 
    
    first=True
    #result  = np.zeros(len(x))
    #
    # Do the basic lines XRF and Escape first...
    #
    for key,linegroup in linedict.iteritems():
        if key in params:
            for line in linegroup:
                area = params[key].value*line.intensity
                center = line.line_energy
                # update funcparams
                if isinstance(line,(XrayLine,EscapeLine)):
                   # for ii in xrf_sigma_keys:
                   #     xrfparams[ii] = np.sqrt(xrfparams[ii]**2 + \
                   #              (center*fano))
                    lineresult =\
                        xrffunc(x=x,area = area,center = center,**xrfparams)
                if isinstance(line,(ElasticLine)):
#                    for ii in elastic_sigma_keys:
#                        params[ii].value = np.sqrt(params[ii].value**2 + \
#                                 (center*fano))
                    lineresult=\
                        elasticfunc(x=x,area = area,center = center,**elasticparams)
                        
                if isinstance(line,(ComptonLine)):
         #           for ii in compton_sigma_keys:
         #               params[ii].value = np.sqrt(params[ii].value**2 + \
         #                        (line.line_energy*0.0001))
                    lineresult=\
                        comptonfunc(x=x,area = area,center = center,**comptonparams)
                if first:
                    result = lineresult
                    first = False
                else:
                    result += lineresult
    #
    # generate the pileup lines...
    # you can convolute the whole line with itself to get the pileup 
    # (basic code is a variation on what's done in pymca..)
    # 
    #
    # The problem is largely in the intensity scaling. For a simple
    # 2-event summation you can just scale it by some factor
    # but for 3 event pileup the magnitude over the curve 
    # grows rather than diminishes.
    #

    #result = result + pileup_factor*selfconvolute(x,result)

    return result
        
        

