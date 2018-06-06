# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model
import numpy as np
from copy import deepcopy



class XRFAllElementModel(Model):

    def __init__(self, exptdesc,linedict,xrffunc,comptonfunc, elasticfunc,
                 linetype='element',*args,**kwargs):
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
        super(XRFAllElementModel, self).__init__(None,*args, **kwargs)
        # Arguments need to be converted to  
        # a list of elements (the peaks -> translating to peak areas) + the line type arguments 
        # The basic model just assigns  
        # self.func = function
        # parses the function parameters to
        #
        # line model
        self.funcmodel = Model(xrffunc)
        self.funcparams_all = self.funcmodel.make_params()
        self.funcparams_sub = self.funcmodel.make_params()   
        # remove the area and center parameter 
        self.funcparams_all.pop("area",None)
        self.funcparams_all.pop("center",None)
        # pileup model
        self.pileupmodel = Model(linefunc)
        self.pileupmodel = self.pileupmodel.make_params()
        # remove the area and center parameter 
        self.funcparams_all.pop("area",None)
        self.funcparams_all.pop("center",None)


        # get the function name
        self.linedict = linedict
        self.exptdesc = exptdesc
        self.funcname   = linefunc.__name__
        self._name = "xrf_model"
        self.func = self.funcover
        #self.create_full_parameters(self)
        
        
    def create_full_parameters(self):
        
        # remove the area and center parameter 
        # For each element in the linedict add a parameter
        # which will control the area
        self.funcparams_all = deepcopy(self.funcparams_sub)
        self.funcparams_all.pop("area",None)
        self.funcparams_all.pop("center",None)
        
        for key in self.exptdesc["lineshape"]["element"][self.funcname].keys():
            if key in self.funcparams_all:
                self.funcparams_all[key].value = \
                    self.exptdesc["lineshape"]["element"]\
                    [self.funcname][key]["value"]
                self.funcparams_all[key].min   = \
                    self.exptdesc["lineshape"]["element"]\
                    [self.funcname][key]["min"]
                self.funcparams_all[key].max   = \
                    self.exptdesc["lineshape"]["element"]\
                    [self.funcname][key]["max"]
                self.funcparams_all[key].vary  = \
                    self.exptdesc["lineshape"]["element"]\
                    [self.funcname][key]["vary"]
        for key in self.linedict:
#            if not '+' in key:            
            self.funcparams_all.add(key,
                **self.exptdesc["lineshape"]["element"][self.funcname]["area"])

         
        

    def make_params(self):
        self.create_full_parameters()
        return self.funcparams_all

    def funcover(self, params,**kwargs):
        #
        # copy the function based parameters from params to funcparams
        # 
        for key in params:
            if key in self.funcparams_sub:
                self.funcparams_sub[key].value = params[key].value
                self.funcparams_sub[key].min   = params[key].min
                self.funcparams_sub[key].max   = params[key].max
                self.funcparams_sub[key].vary  = params[key].vary
        xrf_sigma_keys = [s for s in self.funcparams_sub if 'sigma' in s]
        elastic_sigma_keys = [s for s in self.elasticparams if 'sigma' in s]
        xrf_sigma_keys = [s for s in self.funcparams_sub if 'sigma' in s]
        
        first=True
        #
        # Do the basic lines XRF and Escape first...
        #
        for key,linegroup in self.linedict.iteritems():
            if key in params:
                for line in linegroup:
                    # update funcparams
                    self.funcparams_sub["area"].value   = params[key].value*\
                                                 line.intensity
                    self.funcparams_sub["center"].value = line.line_energy
                    if isinstance(line,(XrayLine,EscapeLine)):
                        for ii in sigma_keys:
                            self.funcparams_sub[ii].value = np.sqrt(params[ii].value**2 + \
                                     (line.line_energy*0.001))
                        lineresult =\
                            self.funcmodel.eval(self.funcparams_sub,**kwargs)
                    if isinstance(line,(ElasticLine)):
                        for ii in elastic_sigma_keys:
                            self.funcparams_sub[ii].value = np.sqrt(params[ii].value**2 + \
                                     (line.line_energy*0.001))
                        lineresult=\
                            self.elasticmodel.eval(self.funcparams_sub,**kwargs)
                    if isinstance(line,(ComptonLine)):
                        for ii in elastic_sigma_keys:
                            self.funcparams_sub[ii].value = np.sqrt(params[ii].value**2 + \
                                     (line.line_energy*0.001))
                        lineresult=\
                            self.comptonmodel.eval(self.funcparams_sub,**kwargs)
                    if isinstance(line,PileupLine):
                            # full method is to convolute the lines 
                            # split the label to determine the individual lines
                            # generate each line...
                            # roll them both to zero
                            # convolve the lineshapes
                            # roll them back to the energy 
                            pass
               
                
                
                    if first:
                        result = lineresult
                        first = False
                    else:
                        result += lineresult
        
        return result
        
        
    def eval(self, params=None, **kwargs):
        result= self.funcover(params=params, **kwargs)
        return result







