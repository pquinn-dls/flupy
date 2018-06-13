# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model,Parameters,minimize
import numpy as np
from flupy.xray.xraylibrary2 import XrayLine,EscapeLine,ComptonLine,ElasticLine
from copy import copy

class XRFAllElementModel():

    def __init__(self, exptdesc,linedict,xrffunc,comptonfunc, elasticfunc,*args,**kwargs):
        """

        An XRF spectrum model
        LMFIT offers a range of options for coupling parameters and models but in order to be able to do 
        this it needs to manipulate parameters - in particular string expressions and this leads 
        big slowdowns if you've a large number of models or coupled parameters.     
        To speed up the calculation rather than create the XRF spectrum from a sum of Gaussian Models
        we create a single model with many parameters.
        The basic lmfit model works as follows:
        parse the arguments of the input function.
        Build a list of param names from these arguments
        Apply paramhints to these param_names or make parameters from the param names
  
        This class will pass a basic line

        """
        self.linedict = linedict
        self.exptdesc = exptdesc

        #
        # xrf line model
        self.xrfmodel = Model(xrffunc,prefix="xrf_",
                              independent_vars=['x','area','center'])
        self.set_param_hints_from_dict(self.xrfmodel,"element",xrffunc.__name__)
        self.xrfparams = self.xrfmodel.make_params()
        #self.xrfparams.pretty_print()
        # elastic model
        self.elasticmodel = Model(elasticfunc,prefix="elastic_",
                                  independent_vars=['x','area','center'])
        self.set_param_hints_from_dict(self.elasticmodel,"elastic",elasticfunc.__name__)
        self.elasticparams = self.elasticmodel.make_params()
        # compton model
        self.comptonmodel = Model(comptonfunc,prefix="compton_",
                                  independent_vars=['x','area','center'])
        self.set_param_hints_from_dict(self.comptonmodel,"compton",comptonfunc.__name__)
        self.comptonparams = self.comptonmodel.make_params()


    def set_param_hints_from_dict(self,model,linetype,linefuncname):

        avail_hints = self.exptdesc["lineshape"][linetype][linefuncname].keys()
        for key in avail_hints:
            if key in model._param_root_names:
                model.set_param_hint(key,expr=None,**self.exptdesc["lineshape"]\
                                      [linetype][linefuncname][key])
        
        
    def make_params(self):
        # A composite of the parameters..
        self.funcparams_all = Parameters()
        self.funcparams_all.update(self.xrfparams)
        self.funcparams_all.update(self.comptonparams)
        self.funcparams_all.update(self.elasticparams)
        
        # exclude pileups from the parameters
        for key in self.linedict:
            if not '+' in key:  
                if "Elastic" in key:
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["elastic"]\
                        [self.elasticmodel.func.__name__]["area"])
                elif "Compton" in key:
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["compton"]\
                        [self.comptonmodel.func.__name__]["area"])
                else:
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["element"]\
                        [self.xrfmodel.func.__name__]["area"])
    #    self.funcparams_all.add("pileup_factor",value=0.05,min=1.0e-5,
    #                           max=2.,vary=True)
                
        return self.funcparams_all
        

    def compositefunc(self, params,**kwargs):

        xrf_sigma_keys     = [s for s in self.xrfparams if 'sigma' in s]
        elastic_sigma_keys = [s for s in self.elasticparams if 'sigma' in s]
        compton_sigma_keys = [s for s in self.comptonparams if 'sigma' in s]
        #
        # common parameters
        # pileup_factor , fano
        #
  #      pileup_factor = params["pileup_factor"].value 
        
        first=True
        x       = kwargs.get('x', None)
        result  = np.zeros(len(x))
        #
        # Do the basic lines XRF and Escape first...
        #
        for key,linegroup in self.linedict.iteritems():
            if key in params:
                for line in linegroup:
                    area = params[key].value*line.intensity
                    center = line.line_energy
                    # update funcparams
                    if isinstance(line,(XrayLine,EscapeLine)):
           #             for ii in xrf_sigma_keys:
           #                 params[ii].value = np.sqrt(params[ii].value**2 + \
           #                          (line.line_energy*0.0001))
                        lineresult =\
                            self.xrfmodel.eval(params=params,area = area, 
                                               center = center,**kwargs)
                    if isinstance(line,(ElasticLine)):
            #            for ii in elastic_sigma_keys:
            #                params[ii].value = np.sqrt(params[ii].value**2 + \
            #                         (line.line_energy*0.0001))
                        lineresult=\
                            self.elasticmodel.eval(params=params,area = area, 
                                                   center = center,**kwargs)
                            
                    if isinstance(line,(ComptonLine)):
             #           for ii in compton_sigma_keys:
             #               params[ii].value = np.sqrt(params[ii].value**2 + \
             #                        (line.line_energy*0.0001))
                        lineresult=\
                            self.comptonmodel.eval(params=params,area = area, 
                                                   center = center,**kwargs)
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

   #     result = result + pileup_factor*self.selfconvolute(result,**kwargs)
    
        return result
        
        
    def eval(self, params=None, **kwargs):
        result= self.compositefunc(params=params, **kwargs)
        return result


    def selfconvolute(self,spectra,**kwargs):
        x = kwargs.get('x',None)
        start_index = int(x[0]/0.01)
        pileup_window = np.zeros(3*len(spectra))
        result   = np.zeros(len(spectra))
        window   = spectra
        spectrum = copy(spectra)
        for i in range(2):
            pileup  = np.convolve(spectrum,window)
            pileup  = pileup/(window.sum())
            pileup_window[start_index:len(pileup)+start_index] += pileup[:]
            spectrum = copy(pileup_window[0:len(spectra)])
            result = result + spectrum
        return result

    def residuals(self,params,*args,**kwargs):
        # Do this to include errors as weights.
        weights = kwargs.get('weights', None)
        data    = kwargs.get('data', None)
        x       = kwargs.get('x', None)
        model  = self.compositefunc(params,x=x)
        resids = model - data
        return resids/weights

    def fit(self,params,*args,**kwargs):
        return minimize(self.residuals, params,*args,**kwargs)




