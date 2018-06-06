# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model
import numpy as np
from copy import deepcopy
from flupy.flupy.xray.xraylibrary2 import XrayLine,EscapeLine,ComptonLine,ElasticLine


class XRFAllElementModel(Model):

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
        # Pass an empty function to model
        super(XRFAllElementModel, self).__init__(None,*args, **kwargs)
        # Arguments need to be converted to  
        # a list of elements (the peaks -> translating to peak areas) + the line type arguments 
        # The basic model just assigns  
        # self.func = function
        # parses the function parameters to
        #
        # xrf line model
        self.xrfmodel = Model(xrffunc,prefix="xrf_")
        self.funcparams_all = self.xrfmodel.make_params()
        self.xrfparams = self.xrfmodel.make_params()   
        # elastic model
        self.elasticmodel = Model(elasticfunc,prefix="elastic_")
        self.elasticparams = self.elasticmodel.make_params()
        # compton model
        self.comptonmodel = Model(comptonfunc,prefix="compton_")
        self.comptonparams = self.comptonmodel.make_params()

        self.linedict = linedict
        self.exptdesc = exptdesc
        self.func = self.funcover
        
    def create_full_parameters(self):
        
        # remove the area and center parameter 
        # For each element in the linedict add a parameter
        # which will control the area
        self.funcparams_all = deepcopy(self.xrfparams)
        self.funcparams_all.pop("xrf_area",None)
        self.funcparams_all.pop("xrf_center",None)
        if [s for s in self.linedict if "compton" in s.lower()]:
            self.funcparams_all.update(self.comptonparams)
            self.funcparams_all.pop("compton_area",None)
            self.funcparams_all.pop("compton_center",None)
        if [s for s in self.linedict if "elastic" in s.lower()]:    
            self.funcparams_all.update(self.elasticparams)
            self.funcparams_all.pop("elastic_area",None)
            self.funcparams_all.pop("elastic_center",None)
        
        #
        # for each model set the parameters...
        #

        for key in self.funcparams_all:
            if key in self.comptonparams:
                ftype = self.comptonmodel.func.__name__
                pdict =self.exptdesc["lineshape"]["compton"][ftype]
                skey  =self.comptonmodel._strip_prefix(key)
                if skey in pdict:
                    self.funcparams_all[key].value = pdict[skey]["value"]
                    self.funcparams_all[key].min = pdict[skey]["min"]
                    self.funcparams_all[key].max = pdict[skey]["max"]
                    self.funcparams_all[key].vary = pdict[skey]["vary"]
            elif key in self.elasticparams:
                ftype = self.elasticmodel.func.__name__
                pdict =self.exptdesc["lineshape"]["elastic"][ftype]
                skey  =self.elasticmodel._strip_prefix(key)
                if skey in pdict:
                    self.funcparams_all[key].value = pdict[skey]["value"]
                    self.funcparams_all[key].min = pdict[skey]["min"]
                    self.funcparams_all[key].max = pdict[skey]["max"]
                    self.funcparams_all[key].vary = pdict[skey]["vary"]
            elif key in self.xrfparams:
                ftype =self.xrfmodel.func.__name__
                pdict =self.exptdesc["lineshape"]["element"][ftype]
                skey  =self.xrfmodel._strip_prefix(key)
                if skey in pdict:
                    self.funcparams_all[key].value = pdict[skey]["value"]
                    self.funcparams_all[key].min   = pdict[skey]["min"]
                    self.funcparams_all[key].max   = pdict[skey]["max"]
                    self.funcparams_all[key].vary  = pdict[skey]["vary"]
        # create a parameter for each line in the line dict...
        # exclude pileups from the parameters
        for key in self.linedict:
            if not '+' in key:  
                if "Elastic" in key:
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["elastic"]\
                        [self.elasticmodel.func.__name__]["area"])
                elif "Compton" in key:
                    print "here",key,self.comptonmodel.func.__name__
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["compton"]\
                        [self.comptonmodel.func.__name__]["area"])
                else:
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["element"]\
                        [self.xrfmodel.func.__name__]["area"])
                    
        #self._func_allargs = self.funcparams_all.keys()
        # create parameters for common value

   def set_param_hints_from_dict(self):
 
       self.set_param_hint()
       
        
    @property
    def param_names(self):
        """Return parameter names for composite model."""
        return self.funcparams_all.keys()

        



    def make_params(self):
        self.create_full_parameters()
        return self.funcparams_all

    def funcover(self, params,**kwargs):
        #
        # copy the function based parameters from params to funcparams
        # 
        for key in params:
            if key in self.xrfparams:
                self.xrfparams[key].value = params[key].value
                self.xrfparams[key].min   = params[key].min
                self.xrfparams[key].max   = params[key].max
                self.xrfparams[key].vary  = params[key].vary
            if key in self.elasticparams:
                self.elasticparams[key].value = params[key].value
                self.elasticparams[key].min   = params[key].min
                self.elasticparams[key].max   = params[key].max
                self.elasticparams[key].vary  = params[key].vary
            if key in self.comptonparams:
                self.comptonparams[key].value = params[key].value
                self.comptonparams[key].min   = params[key].min
                self.comptonparams[key].max   = params[key].max
                self.comptonparams[key].vary  = params[key].vary
        
       # print self.comptonparams["compton_hi_f_tail"].vary,self.comptonparams["compton_hi_f_tail"].value
                
        xrf_sigma_keys     = [s for s in self.xrfparams if 'sigma' in s]
        elastic_sigma_keys = [s for s in self.elasticparams if 'sigma' in s]
        compton_sigma_keys = [s for s in self.comptonparams if 'sigma' in s]
        self.spectradict={}
        first=True
        #
        # Do the basic lines XRF and Escape first...
        #
        for key,linegroup in self.linedict.iteritems():
            if key in params:
                for line in linegroup:
                    # update funcparams
                    if isinstance(line,(XrayLine,EscapeLine)):
                        self.xrfparams["xrf_area"].value   = params[key].value*\
                                                     line.intensity
                        self.xrfparams["xrf_center"].value = line.line_energy

                        for ii in xrf_sigma_keys:
                            self.xrfparams[ii].value = np.sqrt(params[ii].value**2 + \
                                     (line.line_energy*0.001))
                        lineresult =\
                            self.xrfmodel.eval(self.xrfparams,**kwargs)
                    if isinstance(line,(ElasticLine)):
                        self.elasticparams["elastic_area"].value   = params[key].value*\
                                                     line.intensity
                        self.elasticparams["elastic_center"].value = line.line_energy

                        for ii in elastic_sigma_keys:
                            self.elasticparams[ii].value = np.sqrt(params[ii].value**2 + \
                                     (line.line_energy*0.001))
                        lineresult=\
                            self.elasticmodel.eval(self.elasticparams,**kwargs)
                            
                    if isinstance(line,(ComptonLine)):
                        self.comptonparams["compton_area"].value   = params[key].value*\
                                                     line.intensity
                        self.comptonparams["compton_center"].value = line.line_energy
                        for ii in compton_sigma_keys:
                            self.comptonparams[ii].value = np.sqrt(params[ii].value**2 + \
                                     (line.line_energy*0.001))

                        #self.comptonparams.pretty_print()
                        lineresult=\
                            self.comptonmodel.eval(self.comptonparams,**kwargs)
                    if first:
                        result = lineresult
                        first = False
                    else:
                        result += lineresult
                self.spectradict[key] = result
        #
        # generate the pileup lines...
        #        
#        pileup = np.zeros(3*len(result))
#        xmin = 0
#        offset = 0
#        for i in range(len(result)):
#            pileup[i+xmin-offset:i+len(result)+xmin-i-offset] += 1. * result[i] *result[0:len(result)-i]
#        pileup = pileup[0:len(result)]
#        result = result+ pileup
        
        return result
        
        
    def eval(self, params=None, **kwargs):
        result= self.funcover(params=params, **kwargs)
        return result


    def _parse_params(self):
        self._func_haskeywords = (self.left._func_haskeywords or
                                  self.right._func_haskeywords)
        self._func_allargs = (self.left._func_allargs +
                              self.right._func_allargs)
        self.def_vals = deepcopy(self.right.def_vals)
        self.def_vals.update(self.left.def_vals)
        self.opts = deepcopy(self.right.opts)
        self.opts.update(self.left.opts)

 #   def _reprstring(self, long=False):
 #       return "(%s %s %s)" % (self.left._reprstring(long=long),
#                               self._known_ops.get(self.op, self.op),
#                               self.right._reprstring(long=long))

#    def eval(self, params=None, **kwargs):
#        """TODO: docstring in public method."""
#        return self.op(self.left.eval(params=params, **kwargs),
#                       self.right.eval(params=params, **kwargs))

    def eval_components(self, **kwargs):
        """Return OrderedDict of name, results for each component."""
        out = OrderedDict(self.left.eval_components(**kwargs))
        out.update(self.right.eval_components(**kwargs))
        return out


    @property
    def param_names(self):
        """Return parameter names for composite model."""
        return self.left.param_names + self.right.param_names

    @property
    def components(self):
        """Return components for composite model."""
        return self.left.components + self.right.components


    def _get_state(self):
        return (self.left._get_state(),
                self.right._get_state(), self.op.__name__)


    def _set_state(self, state, funcdefs=None):
        return _buildmodel(state, funcdefs=funcdefs)


    def _make_all_args(self, params=None, **kwargs):
        """Generate **all** function arguments for all functions."""
        out = self.right._make_all_args(params=params, **kwargs)
        out.update(self.left._make_all_args(params=params, **kwargs))
        return out





