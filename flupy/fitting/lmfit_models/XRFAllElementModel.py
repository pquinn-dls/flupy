# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 17:52:53 2018

@author: pq67
"""
from lmfit import Model
from itertools import combinations_with_replacement, product
import flupy.flupy.math.lineshapes as lineshapes
import inspect
import numpy as np
from copy import deepcopy
from flupy.flupy.xray.xraylibtools import get_element_xrf_lines,calculate_transmission_corrections

class XRFAllElementModel(Model):

    def __init__(self, paramdict,elementlist,*args,**kwargs):
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
        super(XRFAllElementModel, self).__init__(None, *args, **kwargs)
        # Arguments need to be converted to  
        # a list of elements (the peaks -> translating to peak areas) + the line type arguments 
        # The basic model just assigns  
        # self.func = function
        # parses the function parameters to
        #
        self.paramdict   = paramdict
        self.orig_elementlist = elementlist
        self.elementlist = elementlist

        self.func = self.get_myfunc("split_pvoigt")
        #self.func = self.get_myfunc("gaussian")
        
        self.extract_parameters()        
        
        self._special_parse_params()
        self._name = "xrf_model"

        
    def extract_parameters(self):
        # pull some commonly used values from the metadata
        self.incident_energy        = self.paramdict["Acquisition_instrument"]["XRF"]["beam_energy"]
        self.detector_type          = self.paramdict["Acquisition_instrument"]["XRF"]["Detector"]["detector_type"]
        self.lowE                   = self.paramdict['XRFFit']["FitParams"]["fitted_energy_range_keV"][0]
        self.highE                  = self.paramdict['XRFFit']["FitParams"]["fitted_energy_range_keV"][1]
        self.include_pileup         = self.paramdict['XRFFit']["FitParams"]["include_pileup"]
        self.include_escape         = self.paramdict['XRFFit']["FitParams"]["include_escape"]
        # Build the XRF peak model 
        #self.energies               = np.linspace(self.lowE,self.highE,2000.0)
        self.energies               = self.generate_energy_axis()
        self.transmission_curve     = calculate_transmission_corrections(self.paramdict,self.energies)
        
        self.element_line_dict      = self._generate_element_lines()
        # 
        el_params = set(self.element_line_dict.keys()) & set(self.paramdict["Elements"].keys())
        
        self.pileup_combinations = combinations_with_replacement(el_params,2)


    def generate_energy_axis(self):
        nchannels         = self.paramdict["Acquisition_instrument"]["XRF"]["Detector"]["no_of_channels"]
        gain              = self.paramdict["Acquisition_instrument"]["XRF"]["Detector"]["gain"]
        offset            = self.paramdict["Acquisition_instrument"]["XRF"]["Detector"]["offset"]
        energies          = offset + np.arange(nchannels)*gain
        return energies
   
    def update(self):
        self.extract_parameters()

    def get_myfunc(self,linefunc):
        # Start with the element function...
        func = getattr(lineshapes,linefunc, None)
        if func == None:
            func = getattr(lineshapes,"gaussian", None)
        return func

    
    
    def _special_parse_params(self):
        """Build parameters from function arguments."""
        # Start with the element function....

        # Start with the element function...
        try:  # PY3
            argspec = inspect.getfullargspec(self.func)
            keywords_ = argspec.varkw
        except AttributeError:  # PY2
            argspec = inspect.getargspec(self.func)
            keywords_ = argspec.keywords

        self.pos_args_el = argspec.args
        kw_args = {}
        if argspec.defaults is not None:
            for val in reversed(argspec.defaults):
                kw_args[pos_args.pop()] = val

        # The centre will be the line position so it's not a parameter
        # The area won't be a single parameter - it'll be an area for each element
        #
        # the arguments are now - x, list of elements (areas or concentrations) + lineshape parameters (with area removed)      
        #print pos_args_el[0], self.elementlist, pos_args_el[1:]
        extra_args=[]
        if self.include_pileup:
            extra_args.append("pileup_factor")
        extra_args.append("fwhm_fanoprime")
        pos_args = [self.pos_args_el[0]]+self.elementlist + extra_args+ self.pos_args_el[1:]
        pos_args.remove("area")
        pos_args.remove("center")
        self.pos_args_el.remove("area")
        self.pos_args_el.remove("center")
        self.pos_args_el.pop(0)
        
        self._param_root_names = None
        allargs = pos_args + list(kw_args.keys())
        if len(allargs) == 0 and keywords_ is not None:
            return

        # default independent_var = 1st argument
        self.independent_vars = [pos_args[0]]
        # default param names: all positional args
        # except independent variables
        self.def_vals = {}
        might_be_param = []
        if self._param_root_names is None:
            self._param_root_names = pos_args[:]
            for key, val in kw_args.items():
                if (not isinstance(val, bool) and
                        isinstance(val, (float, int))):
                    self._param_root_names.append(key)
                    self.def_vals[key] = val
                elif val is None:
                    might_be_param.append(key)
            for p in self.independent_vars:
                if p in self._param_root_names:
                    self._param_root_names.remove(p)

        # set the parameter values...
        for mykey in ["pileup","element","common"]:
            for pname in self._param_root_names:
                try:
                    if pname in self.paramdict["XRFFit"][mykey].keys():
                        parm = self.paramdict["XRFFit"][mykey][pname]
                        varyp=self.paramdict["XRFFit"][mykey][pname]["bound_type"]
                        if varyp =='lohi':
                            self.set_param_hint(pname, value=parm["value"],vary=True,min=parm["min"], max=parm["max"])
                        elif varyp=='free':
                            self.set_param_hint(pname, value=parm["value"],vary=True,min=-np.inf, max=np.inf)
                        else:
                            self.set_param_hint(pname, value=parm["value"],vary=False,min=parm["min"], max=parm["max"])
                        continue
                    elname,linename = self.split_elementname(pname)
                    if elname in self.paramdict["Elements"].keys():
                        self.set_param_hint(pname,value=1.0,vary=True,min=0.0,max=np.inf)
                        continue
                    #if pname == "pileup_fraction"
                    #else:
                    #    self.set_param_hint(pname,value=1.0,vary=True,min=-100.0,max=np.inf)
                except KeyError:
                    print("No parameter data found for :",pname)
        new_opts = {}
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

    def split_elementname(self,elementname):
        el = elementname.split('_')
        if len(el)!=2:
            element = el[0]
            lines = 'All'
        else:# len(el)==2:
            element = el[0]
            lines = el[1]
        #elif len(el) > 2:    
        return element,lines        

    def _generate_element_lines(self):
        """
        
        Construct element model. Elemental line that can be add include elemental
        peak, pileup peak and user peak. User peak is a fake peak, but can be used
        as a quick replacement for pileup peaks, escape peaks or other unknown peaks
        due to detector defects or scattering process.

        Parameters
        ----------
        elemental_line : str
        elemental line, such as 'Fe_K'
        pileup peak, such as 'Si_Ka1-Si_Ka1'
        user peak, such as 'user_peak1', 'user_peak2'
        default_area : float, optional
            value for the initial area of a given element
            default is 1e5, found to be a good value
            
        """
        element_line_dict={}
        element_line_dict["xrf"]={}
        element_line_dict["pileup"]={}
        
        for elementitem in self.elementlist:
            # atomic number of the element
            element,lines = self.split_elementname(elementitem)
                
            Z  = self.paramdict["Elements"][element]['Z']

            # get the lines for this element and shell...
            emission_lines = get_element_xrf_lines(Z,lines,self.incident_energy,self.lowE,self.highE,include_escape=self.include_escape,detector_element=self.detector_type)

            # now adjust the emission lines for transmission correction
            for line in emission_lines["lines"]:
                trans_factor =  np.interp(line[1], self.energies,self.transmission_curve)
                line[2] = line[2]*trans_factor

            # now adjust the emission lines for transmission correction
            for line in emission_lines["escape"]:
                trans_factor =  np.interp(line[1], self.energies,self.transmission_curve)
                line[3] = line[3]*trans_factor
            
            # now normalize so that the area of xrf curve will be 1
            # normalize the crossections to 1
            sum_cs = 0.0
            for  line in emission_lines["lines"]:
                sum_cs += line[2]
 
            # scale the cross_sections  
            for line in emission_lines["lines"]:
                line[2] = line[2]/sum_cs

            # scale the escape peaks by the same norm factor
            for line in emission_lines["escape"]:
                line[3] = line[3]/sum_cs


            # if no emission lines due to excitation energy or bounds on fit then...skip to next element
            if not emission_lines:
               # logger.debug('%s emission line is not activated at this energy %f', element, self.incident_energy)
                continue
            else:
                element_line_dict["xrf"][elementitem]=emission_lines

                #element_line_list.append((element,emission_lines))
        if self.include_pileup:
            pileup_combinations = combinations_with_replacement(element_line_dict["xrf"].keys(),2)
            for pileup_combo in pileup_combinations:
                el1=pileup_combo[0]
                el2=pileup_combo[1]
                pileup_lines = self.calculate_pileup_lines(element_line_dict["xrf"][el1]["lines"],
                                                           element_line_dict["xrf"][el2]["lines"])
                # if no emission lines due to excitation energy or bounds on fit then...skip to next element
                if not pileup_lines:
                   # logger.debug('%s %s pileup is not activated at this energy %f or across fitted range %f to %f keV', 
                    #     pileup_combo[0], pileup_combo[1])
                    continue
                element_line_dict["pileup"][(el1,el2)]=pileup_lines
        return element_line_dict

    def calculate_pileup_lines(self,linelist1,linelist2): 
        """
        
        """
        line_combos = product(linelist1,linelist2)
        pileup_list=[]
        for pair in line_combos:
            lineE = pair[0][1]+pair[1][1]
            cs    = pair[0][2]*pair[1][2]
            if(lineE<self.highE):
                pileup_list.append((pair[0][0]+'_'+pair[1][0],lineE,cs))
        return pileup_list                
 
        

    def create_xrf_line(self,params,element,x,norm_to_one=False):
        """
        create an XRF line for the selected element using the 
        data from the model parameters and parameter dictionary
        
        """
        total=np.zeros_like(x)
        fano = 2.96*params["fwhm_fanoprime"].value
        line_args=[]
        indices = [i for i, s in enumerate(self.pos_args_el) if 'sigma' in s]
   #     sum_total = 1.0
        for argv in self.pos_args_el:
             line_args.append(params[argv].value)
        if self.element_line_dict["xrf"][element]:
            line_dict = self.element_line_dict["xrf"][element]  
            par = params[element]
            area    = par.value
            for line_data in line_dict["lines"]:
                center  = line_data[1]
                ratio_v = line_data[2]
                for ii in indices:
                    argv = self.pos_args_el[ii]
                    line_args[ii] = np.sqrt(params[argv].value**2 + (center*fano))
                total_args = [ratio_v,center] + line_args
                total = total + self.func(x,*total_args)
          #  if (norm_to_one):
          #      sum_total = np.sum(total)
          #      total     = total/sum_total
            total = total*area
            if line_dict["escape"]: 
                for line_data in line_dict["escape"]:
                    center  = line_data[2]
                    ratio_v = line_data[3] 
                    for ii in indices:
                        argv = self.pos_args_el[ii]
                        line_args[ii] = np.sqrt(params[argv].value**2 + (center*fano))
                   # total_args =     [area*ratio_v/sum_total,center] + line_args
                    total_args =     [area*ratio_v,center] + line_args
                    total = total + self.func(x,*total_args)

        return total
        
    def create_pileup_line(self,params,el1,el2,x,norm_to_one=False):
        total=np.zeros_like(x)
        fano = 2.96*params["fwhm_fanoprime"].value
        if "pileup_factor" in params.items():
            pileup_fraction = params["pileup_factor"].value
        else:
            pileup_fraction = 1.0
        line_args=[]
        indices = [i for i, s in enumerate(self.pos_args_el) if 'sigma' in s]
        for argv in self.pos_args_el:
             line_args.append(params[argv].value)
        if (el1,el2) in self.element_line_dict["pileup"]:
            line_data = self.element_line_dict["pileup"][(el1,el2)]
            area    = 0.5*(params[el1].value+params[el2].value)
            for pline in line_data:
                center  = pline[1]
                ratio_v = pline[2]
                for ii in indices:
                    argv = self.pos_args_el[ii]
                    line_args[ii] = np.sqrt(params[argv].value**2 + (center*fano))
                total_args =     [ratio_v,center] + line_args
                total = total + self.func(x,*total_args)
            if (norm_to_one):
                sum_total = np.sum(total)
                total  = total/sum_total
            total = total*area*pileup_fraction
        return total
    
    def create_linear_matrix(self,x,params):
        """
        
        Creates a matrix of lineshapes for use with linear fitting (SVD,NNLS) 
        approaches
        
        Returns a  linear matrix and constraint matrix
        
        """
        # generate the element lines
        selected_elements = []
        pileup_elements = []
        matv = []
        kwargs={}
        kwargs["norm_to_one"]=True
        #new_params = deepcopy(params)
        for name, par in params.items():
            if name in self.paramdict["Elements"].keys():
                # calculate the derivative...
                args = [params,name,x]
                deriv =self._numerical_deriv(self.create_xrf_line,args,kwargs)
#                deriv =self._numerical_deriv(x,self.create_xrf_line,params,name)
                
                matv.append(deriv)
                selected_elements.append(name)
                
        # generate the pileup lines

        if self.include_pileup:
            for pileup_combo in self.pileup_combinations:  
                if pileup_combo in self.element_line_dict:
                    #if(pileup_combo[0]==pileup_combo[1]):
                    # calculate the derivative...
                    args = [x,params,pileup_combo[0],pileup_combo[1]]
                    deriv = self._numerical_deriv(self.create_xrf_line,args,kwargs)
#                   deriv =self._numerical_deriv(x,self.create_pileup_line,params,pileup_combo[0],pileup_combo[1])
                        
                    matv.append(deriv)
                    pileup_elements.append((pileup_combo[0],pileup_combo[1]))
                    #else:
                    #    # calculate the derivative...
                    #    args = [pileup_combo[0],pileup_combo[1],norm_to_one=True]
                    #    deriv =self._numerical_deriv(x,self.create_pileup_line,params,args)
                    #    #deriv =self._numerical_deriv(x,self.create_pileup_line,params,pileup_combo[0],pileup_combo[1])
                    #    
                    #    matv.append(deriv)
                    #    pileup_elements.append((pileup_combo[0],pileup_combo[1]))
                    #    args = [pileup_combo[1],pileup_combo[0],norm_to_one=True]
                    #    deriv =self._numerical_deriv(x,self.create_pileup_line,params,args)
                    #   # deriv =self._numerical_deriv(x,self.create_pileup_line,params,pileup_combo[1],pileup_combo[0])
                    #
                    #    matv.append(deriv)
                    #    pileup_elements.append((pileup_combo[1],pileup_combo[0]))
        # spectra matrix
        nspectra = len(matv)
        matv = np.array(matv)
        
        # now generate the constraints matrix...
        nspectra    = len(matv)
        #nelements   = len(selected_elements)
        constraints = np.eye(nspectra,nspectra)

        return matv, constraints, selected_elements+pileup_elements
        
    def _numerical_deriv(self,func,args,kwargs):
        """
        
        numerical derivative...
        
        """
        h = 1.0e-8
        temp = params[args[2]].value
        y = func(*args,**kwargs)
        
        params[args[2]].value = temp + h 
        y_h = func(*args,**kwargs)
        
        params[args[2]].value = temp
        
        deriv = (y_h - y) / h
        return deriv
        
    def eval(self, params=None, **kwargs):
        """Evaluate the model with supplied parameters and keyword arguments.
        Parameters
        -----------
        params : Parameters, optional
            Parameters to use in Model.
        **kwargs : optional
            Additional keyword arguments to pass to model function.
        Returns
        -------
        numpy.ndarray
            Value of model given the parameters and other arguments.
        Notes
        -----
        1. if `params` is None, the values for all parameters are
        expected to be provided as keyword arguments.  If `params` is
        given, and a keyword argument for a parameter value is also given,
        the keyword argument will be used.
        2. all non-parameter arguments for the model function, **including
        all the independent variables** will need to be passed in using
        keyword arguments.
        """
        
        #
        # params is an ordered dict of paremters
        #
#        total=np.zeros_like(x)
        
        #
        # Convert input X data in channels to energy parameters
        #
        
        #
        #
        # XRF peaks including escape peaks
        #
        #if params == None:
        #    params = self.parmas
       # print "**********************************************"
       # print params.pretty_print()
       # print "**********************************************"
        
        first = True
        for name, par in params.items():
            elname,lines = self.split_elementname(name)
            if elname in self.paramdict["Elements"].keys():
                if (first):
                    total = self.create_xrf_line(params,name,**kwargs)
                    first = False
                else:
                    total = total + self.create_xrf_line(params,name,**kwargs)
        
        # Pileup peaks
        #
        if self.include_pileup:
            for pileup_combo in self.pileup_combinations:  
                   if pileup_combo in self.element_line_dict:
                        total=total+self.create_pileup_line(params,pileup_combo[0],pileup_combo[1],**kwargs)
    
        return total


