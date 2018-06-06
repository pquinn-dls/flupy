# ######################################################################
# Copyright (c) 2014, Brookhaven Science Associates, Brookhaven        #
# National Laboratory. All rights reserved.                            #
#                                                                      #
# @author: Li Li (lili@bnl.gov)                                        #
# created on 09/10/2014                                                #
#                                                                      #
# Original code:                                                       #
# @author: Mirna Lerotic, 2nd Look Consulting                          #
#         http://www.2ndlookconsulting.com/                            #
# Copyright (c) 2013, Stefan Vogt, Argonne National Laboratory         #
# All rights reserved.                                                 #
#                                                                      #
# Redistribution and use in source and binary forms, with or without   #
# modification, are permitted provided that the following conditions   #
# are met:                                                             #
#                                                                      #
# * Redistributions of source code must retain the above copyright     #
#   notice, this list of conditions and the following disclaimer.      #
#                                                                      #
# * Redistributions in binary form must reproduce the above copyright  #
#   notice this list of conditions and the following disclaimer in     #
#   the documentation and/or other materials provided with the         #
#   distribution.                                                      #
#                                                                      #
# * Neither the name of the Brookhaven Science Associates, Brookhaven  #
#   National Laboratory nor the names of its contributors may be used  #
#   to endorse or promote products derived from this software without  #
#   specific prior written permission.                                 #
#                                                                      #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS  #
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT    #
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS    #
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE       #
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,           #
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR   #
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)   #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,  #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OTHERWISE) ARISING   #
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE   #
# POSSIBILITY OF SUCH DAMAGE.                                          #
########################################################################
########################################################################
#
# Hacked and modified version of scikit-beam code
#
# Changed
# Made a per element spectrum rather than per emission line...
# lmfit slows down a lot with increased number of models
# A lot of time spent stripping strings
#
# Changed the XRF line creation algorithm and interaction with xraylib
# Also a number of important lines missing.....added back in
#
# Migrated helper and fitter functions to other modules in the hope of
# improving modularity and organization...a bit subjective buy hey. 
#
# integrated escape peaks into the line calculation
#
# Automated pileup peak creation
# 
# Removed the hard coded parameter names - should be able to change the 
# lineshape without editing the underlying code
#
#
#########################################################################
import copy
import logging

from flupy.algorithms.xrf_lmfit_models.compton_model import ComptonModel
from flupy.algorithms.xrf_lmfit_models.elastic_model import ElasticModel
from flupy.algorithms.xrf_lmfit_models.element_model import ElementModel
from flupy.algorithms.xrf_lmfit_models.pileup_model import PileupModel
from flupy.algorithms.xrf_lmfit_models.model_tools import _set_parameter_hint,_link_model_param_hints,_set_parameter_and_link_hint
from flupy.algorithms.xrf_calculations.transitions_and_shells import *
from itertools import combinations_with_replacement,combinations
import xraylib
from lmfit import Model


xraylib.XRayInit()
logger = logging.getLogger(__name__)

class XRFModelSpectrum(object):
    """

    Construct Fluorescence spectrum which includes elastic peak,
    compton and a single element peaks. To be used for fitting
    a XRF fitting a series of spectra.

    """

    def __init__(self, params, elementlist):
        """
        Parameters
        ----------
        params : dict of lmfit params
            saving all the fitting values and their bounds
        elemental -  e.g. Fe
        shell     -  e.g K or L1 , L2, L3 etc,

        """
        self.params            = params
        # pull some values from the fit configurations
        self.incident_energy        = self.params['XRFFit']["common"]['coherent_sct_energy']['value']
        self.lowE                   = self.params['XRFFit']["FitParams"]["fitted_energy_range_keV"][0]
        self.highE                  = self.params['XRFFit']["FitParams"]["fitted_energy_range_keV"][1]
        self.include_pileup         = self.params['XRFFit']["FitParams"]["include_pileup"]

        self.modelDict = {}
        self.params = params
        self.elementModelList = []
        self.pileupModelList = []

        self.setup_compton_model()
        self.setup_elastic_model()
        # Build the XRF peak model 
        self.elementlist = elementlist
        self.setup_element_model()
        #self.setup_pileup_model()
        if(self.include_pileup):
            self.setup_pileup_model()

        # bring them together and clean up the share of core parameters
        self.assembleFullModel()

    def setup_compton_model(self):
        """
        setup parameters related to Compton model

        """
        # create the model
        self.modelDict["compton"]=ComptonModel(prefix='compton_',name='compton')


    def setup_elastic_model(self):
        """
        setup parameters related to Elastic model
        """

        # create the model
        self.modelDict["elastic"]=ElasticModel(prefix='elastic_',name='elastic')


    def setup_element_model(self,default_area=1e5):
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
        self.modelDict["element"]={}
        for element in self.elementlist:
            # atomic number of the element
            Z  = self.params["Elements"][element]['Z']

            # get the lines for this element and shell...
            emission_lines = xl_emmission_lines_energy(Z,self.incident_energy,self.lowE,self.highE)

            # if no emission lines due to excitation energy or bounds on fit then...skip to next element
            if not emission_lines:
                logger.debug('%s emission line is not activated at this energy %f', element, self.incident_energy)
                continue

            lineprefix=str(element)+'_'
            element_mod = ElementModel(prefix=lineprefix,emission_lines=emission_lines,line_type="element")
            self.modelDict["element"][element]=element_mod
            logger.debug('Finished building element peak for %s', element)

    def assembleFullModel(self):
	#
        # self.modeldict holds the models for elastic,compton,elements and pileups

        # 1) Set the element default parameter values
        # 2) Set the compton default parameter values
        # 3) Set the elastic default parameter values
        # 4) Set the pileup default parameter values

        # 2) link the element peak shape parameters (except area)
        # 5) link energy scaling or common parameters of all models
        # 5) link scattering energy common parameters of all models
        #    link pileup areas to element areas
        
        # 1-4 set the defaults
        self.set_defaults_from_config(self.modelDict)      

        # 5- set the element links
        current_element_model=self.modelDict["element"].itervalues().next()
        base_prefix= current_element_model.prefix
        # link the element parameters
        self.set_links_from_config(self.modelDict,expr1="",expr2=base_prefix,models_to_link=["element"])
        current_element_model=self.modelDict["pileup"].itervalues().next()
        base_prefix= current_element_model.prefix
        self.set_links_from_config(self.modelDict,expr1="",expr2=base_prefix,models_to_link=["pileup"])
        # 6-7 link the common parameters to the compton model      
        self.set_links_from_config(self.modelDict,expr1="",expr2="compton_",models_to_link=["all"]) 
        # 8 link the pileup area parameters to the element model      

        self.mod = self.sumModels(self.modelDict)
        # add the pileup fraction
        self.all_params = self.mod.make_params()
        self.all_params.pretty_print()


    def setup_pileup_model(self,default_area=1e5):
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
        self.modelDict["pileup"]={}
        # pull some values from the fit configurations
        epsilon = self.params["XRFFit"]["non_fitting_values"]["epsilon"]

        all_pileup_mod = None
        # global parameters - smae energy axis and fwhm parameters.
        firstelement=False

        # looking at 2 event pileup
        # get the possible combinations of elements 
        pileup_combinations = combinations_with_replacement(self.elementlist,2)
        for pileup_combo in pileup_combinations:
            # atomic number of the element
            Z1  = self.params["Elements"][pileup_combo[0]]['Z']
            Z2  = self.params["Elements"][pileup_combo[1]]['Z']

            # get the lines for this element and shell...
            pileup_lines = xl_emmission_pileup_lines_energy(Z1,Z2,self.incident_energy,self.lowE,self.highE)

            # if no emission lines due to excitation energy or bounds on fit then...skip to next element
            if not pileup_lines:
                logger.debug('%s %s pileup is not activated at this energy %f or across fitted range %f to %f keV', 
                             pileup_combo[0], pileup_combo[1],self.incident_energy,self.lowE,self.highE)
                continue
            
            line_name = pileup_lines[0][1]

            prefix=str(pileup_combo[0])+'_'+str(pileup_combo[1])+"_"
            firstelement_prefix=prefix
            pileup_mod = PileupModel(prefix=prefix,name="pileup",emission_lines=pileup_lines)
            self.modelDict["pileup"][prefix]=pileup_mod
 
    def set_defaults_from_config(self,mydict):
      #
      # recursive loop to set the defaults for all models  
      #
      for k,v in mydict.iteritems():
          if isinstance(v, dict):
              self.set_defaults_from_config(v)
          else:
              if isinstance(v, ElementModel):
                  model_type="element"
              elif isinstance(v, ComptonModel):
                  model_type="compton"
              elif isinstance(v,ElasticModel):
                  model_type="elastic"
              elif isinstance(v,PileupModel):
                  model_type="pileup"
              for item in v.param_names:
                # get the param name and strip the prefix
                sitem = self.param_name_strip(v.prefix,item)
                # check if we've defined defaults for this...
                # set and link model specific parameters
                if sitem in self.params["XRFFit"][model_type].keys():
                    _set_parameter_hint(sitem, self.params["XRFFit"][model_type][sitem], v)
                if sitem in self.params["XRFFit"]["common"].keys():
                    _set_parameter_hint(sitem, self.params["XRFFit"]["common"][sitem], v)

    def sumModels(self,mydict):
        total = None
        for key, value in mydict.iteritems():
            if isinstance(value, Model):
                if total:
                    total += value
                else:
                    total = value
            if isinstance(value, dict):
                if total: 
                    total += self.sumModels(value)
                else:
                    total = self.sumModels(value)
        return total

    def set_links_from_config(self,mydict,expr1="",expr2="",models_to_link=["elements","compton","pileup","elastic"]):
      #
      # recursive loo[ to set the links between modesl
      # link = expr1+expr2+parameter_name
      # normally expr1="" , expr= prefix of paramter to link to
      # 
       for k,v in mydict.iteritems():
          if isinstance(v, dict):
              self.set_links_from_config(v,expr1=expr1,expr2=expr2,models_to_link=models_to_link)
          else:
              # if the baseprefix then this is the one the others are linked to...
              if(v.prefix==expr2 or v.prefix==expr1):
                 continue
              if isinstance(v, ElementModel):
                  model_type="element"
              elif isinstance(v, ComptonModel):
                  model_type="compton"
              elif isinstance(v,ElasticModel):
                  model_type="elastic"
              elif isinstance(v,PileupModel):
                  model_type="pileup"
              for item in v.param_names:
                sitem = self.param_name_strip(v.prefix,item)
                # don't link areas
                # case 1 - fwhm of all elements linked but not lnked to fwhm of other line types..
                if(sitem=="area" and model_type in ["element","elastic","compton"]):
                   continue
                if (model_type in models_to_link):
                    # special cases is pileup where the areas need to be linked
                    # set and link model specific parameters
                    if(model_type=="pileup" and sitem=="area"):
                        els = v.prefix.split("_")  
                        t1=els[0]+"_area+"
                        t2=els[1]+"_area"
                        _set_parameter_and_link_hint(sitem,self.params["XRFFit"][model_type][sitem],t1+t2,v)
                    else:
                        if sitem in self.params["XRFFit"][model_type].keys():
                            _set_parameter_and_link_hint(sitem,self.params["XRFFit"][model_type][sitem],expr1+expr2+sitem,v)
                if "all" in models_to_link:
                    if sitem in self.params["XRFFit"]["common"].keys():
                        _set_parameter_and_link_hint(sitem,self.params["XRFFit"]["common"][sitem],expr1+expr2+sitem,v)

    def assemble_models(self):
        """
        Put all models together to form a spectrum.
        """
            
        self.link_models_and_parameters()


    def param_name_strip(self,prefix,name):
        """
            strips the prefix off...as per Model class
        """ 
        npref = len(prefix)
        if npref > 0 and name.startswith(prefix):
            name = name[npref:] 
        return name



    def rebuild_models(self):
        self.incident_energy        = self.params['XRFFit']['coherent_sct_energy']['value']
        self.lowE                   = self.params['XRFFit']["fitted_energy_range_keV"][0]
        self.highE                  = self.params['XRFFit']["fitted_energy_range_keV"][1]
        self.elastic_compton_cutoff = self.params['XRFFit']["elastic_compton_cutoff"]
        self.include_pileup         = self.params['XRFFit']["include_pileup"]
        self.setup_compton_model()
        self.setup_elastic_model()
        self.setup_element_model()
        self.assemble_models()


    def model_fit(self, channel_number, spectrum, weights=None,
                  method='leastsq', **kwargs):
        """
        Parameters
        ----------
        channel_number : array
            independent variable
        spectrum : array
            intensity
        weights : array, optional
            weight for fitting
        method : str
            default as leastsq
        kwargs : dict
            fitting criteria, such as max number of iteration

        Returns
        -------
        result object from lmfit
        """

        pars = self.mod.make_params()
        result = self.mod.fit(spectrum, self.all_params, x=channel_number, weights=weights,
                              method=method, fit_kws=kwargs)
        return result






