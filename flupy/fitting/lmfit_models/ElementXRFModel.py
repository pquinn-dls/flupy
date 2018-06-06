# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model,Parameters
import numpy as np
from flupy.flupy.xray.xraylibrary2 import XrayLine,EscapeLine
from math import sqrt
    
class XRFModel(Model):

    def __init__(self, exptdesc,linedict,linefunc,*args,**kwargs):
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
        # Pass an empty function to model
        super(XRFModel, self).__init__(linefunc,prefix="xrf_",
                              independent_vars=['x','area','center'])
        
        self.linedict = linedict
        self.exptdesc = exptdesc
        self.func = self.xrffunc
        self.linefunc = linefunc
        self.xrfmodel = Model(linefunc)
        self.set_param_hints_from_dict("element",linefunc.__name__)

        name = None
        if name is None and hasattr(self.xrffunc, '__name__'):
            name = self.func.__name__
        self._name = name
        


    def set_param_hints_from_dict(self,linetype,linefuncname):
        """

        This sets the default paramter hints for the model on setup
        
        """
        avail_hints = self.exptdesc["lineshape"][linetype][linefuncname].keys()
        for key in avail_hints:
            if key in self._param_root_names:
                self.set_param_hint(key,expr=None,**self.exptdesc["lineshape"]\
                                      [linetype][linefuncname][key])
        # exclude pileups from the parameters
        for key in self.linedict:
            if not any(word in key for word in ['+','Elastic','Compton']):
                self.set_param_hint(key,expr=None,**self.exptdesc["lineshape"]["element"]\
                    [self.linefunc.__name__]["area"])

        self.set_param_hint("fano",expr=None,**self.exptdesc["xrf_fitting"]["params"]\
                    ["fano"])


    def make_params(self):
        
        params = super(XRFModel, self).make_params()
        #
        # If the element list is empty return no parameters...
        #
        elementlist = []
        for key in self.linedict:
            if any(keyv not in key for keyv in ['+','Elastic','Compton']):
                elementlist.append(key)                    
        if elementlist:
            return params
        else:
            return Parameters()

                
    def xrffunc(self, params,**kwargs):
         
        xrf_sigma_keys     = [s for s in params if 'sigma' in s]
        #
        # common parameters
        # pileup_factor , fano
        #
        first   = True
        x       = kwargs.get('x', None)
        result  = np.zeros(len(x))
        fano = params["fano"]
        #
        # Do the basic lines XRF and Escape first...
        #
        for key,linegroup in self.linedict.iteritems():
            if key in params:
                for line in linegroup:
                    area = params[key].value*line.intensity
                    center = line.line_energy
                    # update funcparams
                    #cif type(line) == XrayLine or type(line) == EscapeLine:
                    if isinstance(line,(XrayLine,EscapeLine)):
                        for ii in xrf_sigma_keys:
                            params[ii].value = sqrt(params[ii].value**2 + \
                                     (center*fano))
                        lineresult =\
                            self.xrfmodel.eval(params=params,area = area, 
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

#        result = result + pileup_factor*self.selfconvolute(result,**kwargs)
    
        return result

    def eval(self, params=None, **kwargs):
        pardict = self.make_funcargs(params)
        return self.xrffunc(params,pardict, **kwargs)


