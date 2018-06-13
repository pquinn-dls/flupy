# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model,Parameters,minimize
import numpy as np
from flupy.xray.xraylibrary2 import XrayLine,EscapeLine,ComptonLine,ElasticLine
from copy import copy

class XRFAllElementModel(Model):

    def __init__(self, exptdesc,linedict,xrffunc,comptonfunc,elasticfunc,*args,**kwargs):
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
        # a special version of parse params..   
        self._parse_params()

    def set_param_hints_from_dict(self,model,linetype,linefuncname):
        """

        This sets the default paramter hints for the model on setup
        
        """
        avail_hints = self.exptdesc["lineshape"][linetype][linefuncname].keys()
        for key in avail_hints:
            if key in model._param_root_names:
                model.set_param_hint(key,expr=None,**self.exptdesc["lineshape"]\
                                      [linetype][linefuncname][key])

    def set_additional_hints(self):
        
        # if xrf elements are present add fano parameter
        
        # adjust energy gain and offset are channels present
        
        # if include pileup add pileup factor as a parameter
        self.funcparams_all.add("pileup_factor",value=0.05,min=1.0e-5,
                               max=2.,vary=True)
        
        
    def make_params(self):
        """Create a Parameters object for this XRF Model.

        Returns
        ---------
        params : Parameters

        Notes
        -----
        1. The parameters may or may not have decent initial values for each
        parameter.

        2. This applies any default values or parameter hints that may have
        been set.

        """
        
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
                    print "here",key,self.comptonmodel.func.__name__
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["compton"]\
                        [self.comptonmodel.func.__name__]["area"])
                else:
                    self.funcparams_all.add(key,
                        **self.exptdesc["lineshape"]["element"]\
                        [self.xrfmodel.func.__name__]["area"])
        self.funcparams_all.add("pileup_factor",value=0.05,min=1.0e-5,
                               max=2.,vary=True)
                         
        return self.funcparams_all
        

    def compositefunc(self, params,**kwargs):
         
        xrf_sigma_keys     = [s for s in self.xrfparams if 'sigma' in s]
        elastic_sigma_keys = [s for s in self.elasticparams if 'sigma' in s]
        compton_sigma_keys = [s for s in self.comptonparams if 'sigma' in s]
        #
        # common parameters
        # pileup_factor , fano
        #
        pileup_factor = params["pileup_factor"].value 
        
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

        result = result + pileup_factor*self.selfconvolute(result,**kwargs)
    
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

    def _parse_params(self):
        self._func_haskeywords = (self.left._func_haskeywords or
                                  self.right._func_haskeywords)
        self._func_allargs = (self.left._func_allargs +
                              self.right._func_allargs)
        self.def_vals = deepcopy(self.right.def_vals)
        self.def_vals.update(self.left.def_vals)
        self.opts = deepcopy(self.right.opts)
        self.opts.update(self.left.opts)
        
        
   def _parse_params(self):
       
        """Build parameters from function arguments."""
        
        if self.func is None:
            return
        if hasattr(self.func, 'argnames') and hasattr(self.func, 'kwargs'):
            pos_args = self.func.argnames[:]
            kw_args = {}
            for name, defval in self.func.kwargs:
                kw_args[name] = defval
            keywords_ = list(kw_args.keys())
        else:
            try:  # PY3
                argspec = inspect.getfullargspec(self.func)
                keywords_ = argspec.varkw
            except AttributeError:  # PY2
                argspec = inspect.getargspec(self.func)
                keywords_ = argspec.keywords

            pos_args = argspec.args
            kw_args = {}
            if argspec.defaults is not None:
                for val in reversed(argspec.defaults):
                    kw_args[pos_args.pop()] = val

        self._func_haskeywords = keywords_ is not None
        self._func_allargs = pos_args + list(kw_args.keys())
        allargs = self._func_allargs

        # default independent_var = 1st argument
        if self.independent_vars is None:
            self.independent_vars = [pos_args[0]]

        # default param names: all positional args
        # except independent variables
        self.def_vals = {}
        might_be_param = []
        if self._param_root_names is None:
            self._param_root_names = pos_args[:]
            
                    self._param_root_names.append(key)
                    self.def_vals[key] = val

            for p in self.independent_vars:
                if p in self._param_root_names:
                    self._param_root_names.remove(p)

        new_opts = {}
        for opt, val in self.opts.items():
            if (opt in self._param_root_names or opt in might_be_param and
                    isinstance(val, Parameter)):
                self.set_param_hint(opt, value=val.value,
                                    min=val.min, max=val.max, expr=val.expr)
            elif opt in self._func_allargs:
                new_opts[opt] = val
        self.opts = new_opts

        names = []
        if self._prefix is None:
            self._prefix = ''
        for pname in self._param_root_names:
            names.append("%s%s" % (self._prefix, pname))

        # check variables names for validity
        # The implicit magic in fit() requires us to disallow some
        fname = self.func.__name__
        for arg in self.independent_vars:
            if arg not in allargs or arg in self._forbidden_args:
                raise ValueError(self._invalid_ivar % (arg, fname))
        for arg in names:
            if (self._strip_prefix(arg) not in allargs or
                    arg in self._forbidden_args):
                raise ValueError(self._invalid_par % (arg, fname))
        # the following as been changed from OrderedSet for the time being.
        self._param_names = names[:]



    def _reprstring(self, long=False):
        return "(%s %s %s)" % (self.left._reprstring(long=long),
                               self._known_ops.get(self.op, self.op),
                               self.right._reprstring(long=long))

    def eval(self, params=None, **kwargs):
        """TODO: docstring in public method."""
        return self.op(self.left.eval(params=params, **kwargs),
                       self.right.eval(params=params, **kwargs))

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

