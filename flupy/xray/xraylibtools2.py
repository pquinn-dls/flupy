# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 09:30:36 2018

@author: pq67
"""

#from __future__ import absolute_import, division, print_function
from collections import namedtuple
from flupy.flupy.xray.escape import calc_escape_peak_ratios,calc_escape_peak_energy
import xraylib
from xraylib import *
import numpy as np
import re
import math
import logging

xraylib.XRayInit()
xraylib.SetErrorMessages(0)
    
logger = logging.getLogger(__name__)

#
# useful containers for line information
#

XrayLine = namedtuple('XrayLine', ('label','line_energy', 'cs'))
EscapeLine = namedtuple('EscapeLine', ('label','line_energy','escape_energy', 'cs'))
PileupLine = namedtuple('PileupLine', ('label','xraylibline','energy', 'cs'))

logger = logging.getLogger(__name__)

siegbahn_name = ['Ka1', 'Ka2', 'Kb1', 'Kb2', 'La1', 'La2', 'Lb1', 'Lb2',
             'Lb3', 'Lb4', 'Lb5', 'Lg1', 'Lg2', 'Lg3', 'Lg4', 'Ll',
             'Ln', 'Ma1', 'Ma2', 'Mb', 'Mg']

bindingE = ['K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1',
            'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3',
            'O4', 'O5', 'P1', 'P2', 'P3']

K_TRANSITIONS = ['ka1', 'ka2', 'kb1', 'kb2','kb3','kb4','kb5']

L_TRANSITIONS = ['la1', 'la2', 'lb1', 'lb2', 'lb3', 'lb4', 'lb5',
                 'lg1', 'lg2', 'lg3', 'lg4', 'll', 'ln']

M_TRANSITIONS = ['ma1', 'ma2', 'mb', 'mg']

TRANSITIONS_LOOKUP = {'K': K_TRANSITIONS, 'L': L_TRANSITIONS,
                      'M': M_TRANSITIONS}

siegbahn_all_list = [KA1_LINE, KA2_LINE, KA3_LINE,
        KB1_LINE,KB2_LINE,KB3_LINE,
        KB4_LINE, KB5_LINE,LA1_LINE, LA2_LINE,
        LB1_LINE, LB2_LINE, LB3_LINE,
        LB4_LINE, LB5_LINE, LB6_LINE,LB7_LINE,LB9_LINE,LB10_LINE,LB15_LINE,LB17_LINE,
        LG1_LINE,LG2_LINE, LG3_LINE,LG4_LINE,LG5_LINE,LG6_LINE,LG8_LINE,
        LL_LINE, LE_LINE, MA1_LINE,
        MA2_LINE, MB_LINE, MG_LINE]

siegbahn_all_name_list=['ka1', 'ka2', 'ka3','kb1', 'kb2', 'kb3','kb4','kb5',
                        'la1', 'la2', 'lb1', 'lb2', 'lb3', 'lb4', 'lb5',
                 'lg1', 'lg2', 'lg3', 'lg4', 'll', 'ln','ma1','ma2','mb','mg']

siegbahn_ka_list = [KA1_LINE, KA2_LINE, KA3_LINE]
siegbahn_ka_name_list = ['ka1', 'ka2', 'ka3']

siegbahn_kb_list = [KB1_LINE,KB2_LINE,KB3_LINE,KB4_LINE, KB5_LINE]
siegbahn_kb_name_list = ['kb1', 'kb2', 'kb3','kb4','kb5']

siegbahn_k_list = [KA1_LINE, KA2_LINE, KA3_LINE,
        KB1_LINE,KB2_LINE,KB3_LINE,
        KB4_LINE, KB5_LINE]
siegbahn_k_name_list = ['ka1', 'ka2', 'ka3','kb1', 'kb2', 'kb3','kb4','kb5']


siegbahn_l_list = [LA1_LINE, LA2_LINE,
        LB1_LINE, LB2_LINE, LB3_LINE,
        LB4_LINE, LB5_LINE, LB6_LINE,
        LG1_LINE, LG2_LINE, LG3_LINE, LG6_LINE,
        LG8_LINE,LL_LINE, LE_LINE,LL_LINE,]
siegbahn_l_name_list = ['la1', 'la2', 'lb1', 'lb2', 'lb3', 'lb4', 'lb5',
                 'lg1', 'lg2', 'lg3', 'lg4', 'll', 'ln']

siegbahn_l1_list = [LA1_LINE, LA2_LINE]
siegbahn_l1_name_list = ['la1', 'la2']

siegbahn_l2_list = [LB1_LINE, LB2_LINE, LB3_LINE,LB4_LINE, LB5_LINE, LB6_LINE,LB7_LINE,LB9_LINE,LB10_LINE,LB15_LINE,LB17_LINE]
siegbahn_l2_name_list =['lb1', 'lb2', 'lb3', 'lb4', 'lb5', 'lb6', 'lb7', 'lb9', 'lb10', 'lb15', 'lb17']

siegbahn_l3_list = [LG1_LINE,LG2_LINE, LG3_LINE,LG4_LINE,LG5_LINE,LG6_LINE,LG8_LINE]
siegbahn_l3_name_list =['lg1', 'lg2', 'lg3', 'lg4', 'lg5', 'lg6', 'lg8']

siegbahn_m_list = [MA1_LINE, MA2_LINE, MB_LINE, MG_LINE]
siegbahn_m_name_list = ['ma1','ma2','mb','mg']

siegbahn_ma_list = [MA1_LINE, MA2_LINE]
siegbahn_ma_name_list = ['ma1','ma2']

siegbahn_mb_list = [MB_LINE]
siegbahn_mb_name_list = ['mb']

siegbahn_mg_list = [MG_LINE]
siegbahn_mg_name_list = ['mg']

#define MA1_LINE M5N7_LINE
#define MA2_LINE M5N6_LINE
#define MB_LINE M4N6_LINE
#define MG_LINE M3N5_LINE

shell_list = [K_SHELL,  L1_SHELL, L2_SHELL,
        L3_SHELL, M1_SHELL, M2_SHELL,
        M3_SHELL, M4_SHELL, M5_SHELL,
        N1_SHELL, N2_SHELL, N3_SHELL,
        N4_SHELL, N5_SHELL, N6_SHELL,
        N7_SHELL, O1_SHELL, O2_SHELL,
        O3_SHELL, O4_SHELL, O5_SHELL,
        P1_SHELL, P2_SHELL, P3_SHELL]

iupac_list=[KL1_LINE,KL2_LINE,KL3_LINE,KM1_LINE,KM2_LINE,
KM3_LINE,KM4_LINE,KM5_LINE,KN1_LINE, KN2_LINE, 
KN3_LINE,KN4_LINE, KN5_LINE,KN6_LINE,KN7_LINE,
KO_LINE,KO1_LINE,KO2_LINE,KO3_LINE,KO4_LINE, 
KO5_LINE,KO6_LINE,KO7_LINE,KP_LINE,
KP1_LINE,KP2_LINE,KP3_LINE,KP4_LINE, 
KP5_LINE,L1L2_LINE,L1L3_LINE,L1M1_LINE, 
L1M2_LINE,L1M3_LINE,L1M4_LINE,L1M5_LINE, 
L1N1_LINE,L1N2_LINE,L1N3_LINE,L1N4_LINE,
L1N5_LINE,L1N6_LINE,L1N67_LINE,L1N7_LINE, 
L1O1_LINE,L1O2_LINE,L1O3_LINE,L1O4_LINE,
L1O45_LINE,L1O5_LINE,L1O6_LINE,L1O7_LINE,
L1P1_LINE,L1P2_LINE,L1P23_LINE,L1P3_LINE,
L1P4_LINE,L1P5_LINE,L2L3_LINE,L2M1_LINE,
L2M2_LINE,L2M3_LINE,L2M4_LINE,L2M5_LINE,
L2N1_LINE,L2N2_LINE,L2N3_LINE,L2N4_LINE,
L2N5_LINE,L2N6_LINE,L2N7_LINE,L2O1_LINE, 
L2O2_LINE,L2O3_LINE, L2O4_LINE, L2O5_LINE, 
L2O6_LINE, L2O7_LINE, L2P1_LINE, L2P2_LINE,
L2P23_LINE,L2P3_LINE,L2P4_LINE,L2P5_LINE,
L2Q1_LINE,L3M1_LINE,L3M2_LINE,L3M3_LINE,
L3M4_LINE,L3M5_LINE,L3N1_LINE,L3N2_LINE,
L3N3_LINE,L3N4_LINE,L3N5_LINE,L3N6_LINE,
L3N7_LINE,L3O1_LINE,L3O2_LINE,L3O3_LINE,
L3O4_LINE,L3O45_LINE,L3O5_LINE,L3O6_LINE,
L3O7_LINE,L3P1_LINE,L3P2_LINE,L3P23_LINE,
L3P3_LINE,L3P4_LINE,L3P45_LINE,L3P5_LINE,
L3Q1_LINE,M1M2_LINE,M1M3_LINE,M1M4_LINE,
M1M5_LINE,M1N1_LINE,M1N2_LINE,M1N3_LINE,
M1N4_LINE,M1N5_LINE,M1N6_LINE,M1N7_LINE,
M1O1_LINE,M1O2_LINE,M1O3_LINE,M1O4_LINE,
M1O5_LINE,M1O6_LINE,M1O7_LINE,M1P1_LINE,
M1P2_LINE,M1P3_LINE,M1P4_LINE,M1P5_LINE,
M2M3_LINE,M2M4_LINE,M2M5_LINE,M2N1_LINE,
M2N2_LINE,M2N3_LINE,M2N4_LINE,M2N5_LINE,
M2N6_LINE,M2N7_LINE,M2O1_LINE,M2O2_LINE,
M2O3_LINE,M2O4_LINE,M2O5_LINE,M2O6_LINE,
M2O7_LINE,M2P1_LINE,M2P2_LINE,M2P3_LINE,
M2P4_LINE,M2P5_LINE,M3M4_LINE,M3M5_LINE,
M3N1_LINE,M3N2_LINE,M3N3_LINE,M3N4_LINE,
M3N5_LINE,M3N6_LINE,M3N7_LINE,M3O1_LINE,
M3O2_LINE,M3O3_LINE,M3O4_LINE,M3O5_LINE,
M3O6_LINE,M3O7_LINE,M3P1_LINE,M3P2_LINE,
M3P3_LINE,M3P4_LINE,M3P5_LINE,M3Q1_LINE,
M4M5_LINE,M4N1_LINE,M4N2_LINE,M4N3_LINE,
M4N4_LINE,M4N5_LINE,M4N6_LINE,M4N7_LINE,
M4O1_LINE,M4O2_LINE,M4O3_LINE,M4O4_LINE,
M4O5_LINE,M4O6_LINE,M4O7_LINE,M4P1_LINE,
M4P2_LINE,M4P3_LINE,M4P4_LINE,M4P5_LINE,
M5N1_LINE,M5N2_LINE,M5N3_LINE,M5N4_LINE,
M5N5_LINE,M5N6_LINE,M5N7_LINE,M5O1_LINE,
M5O2_LINE,M5O3_LINE,M5O4_LINE,M5O5_LINE,
M5O6_LINE,M5O7_LINE,M5P1_LINE,M5P2_LINE,
M5P3_LINE,M5P4_LINE,M5P5_LINE,N1N2_LINE,
N1N3_LINE,N1N4_LINE,N1N5_LINE,N1N6_LINE,
N1N7_LINE,N1O1_LINE,N1O2_LINE,N1O3_LINE,
N1O4_LINE,N1O5_LINE,N1O6_LINE,N1O7_LINE,
N1P1_LINE,N1P2_LINE,N1P3_LINE,N1P4_LINE,
N1P5_LINE,N2N3_LINE,N2N4_LINE,N2N5_LINE,
N2N6_LINE,N2N7_LINE,N2O1_LINE,N2O2_LINE,
N2O3_LINE,N2O4_LINE,N2O5_LINE,N2O6_LINE,
N2O7_LINE,N2P1_LINE,N2P2_LINE,N2P3_LINE,
N2P4_LINE,N2P5_LINE,N3N4_LINE,N3N5_LINE,
N3N6_LINE,N3N7_LINE,N3O1_LINE,N3O2_LINE,
N3O3_LINE,N3O4_LINE,N3O5_LINE,N3O6_LINE,
N3O7_LINE,N3P1_LINE,N3P2_LINE,N3P3_LINE,
N3P4_LINE,N3P5_LINE,N4N5_LINE,N4N6_LINE,
N4N7_LINE,N4O1_LINE,N4O2_LINE,N4O3_LINE,
N4O4_LINE,N4O5_LINE,N4O6_LINE,N4O7_LINE,
N4P1_LINE,N4P2_LINE,N4P3_LINE,N4P4_LINE,
N4P5_LINE,N5N6_LINE,N5N7_LINE,N5O1_LINE,
N5O2_LINE,N5O3_LINE,N5O4_LINE,N5O5_LINE,
N5O6_LINE,N5O7_LINE,N5P1_LINE,N5P2_LINE,
N5P3_LINE,N5P4_LINE,N5P5_LINE,N6N7_LINE,
N6O1_LINE,N6O2_LINE,N6O3_LINE,N6O4_LINE,
N6O5_LINE,N6O6_LINE,N6O7_LINE,N6P1_LINE,
N6P2_LINE,N6P3_LINE,N6P4_LINE,N6P5_LINE,
N7O1_LINE,N7O2_LINE,N7O3_LINE,N7O4_LINE,
N7O5_LINE,N7O6_LINE,N7O7_LINE,N7P1_LINE,
N7P2_LINE,N7P3_LINE,N7P4_LINE,N7P5_LINE,
O1O2_LINE,O1O3_LINE,O1O4_LINE,O1O5_LINE,
O1O6_LINE,O1O7_LINE,O1P1_LINE,O1P2_LINE,
O1P3_LINE,O1P4_LINE,O1P5_LINE,O2O3_LINE,
O2O4_LINE,O2O5_LINE,O2O6_LINE,O2O7_LINE,
O2P1_LINE,O2P2_LINE,O2P3_LINE,O2P4_LINE,
O2P5_LINE,O3O4_LINE,O3O5_LINE,O3O6_LINE,
O3O7_LINE,O3P1_LINE,O3P2_LINE,O3P3_LINE,
O3P4_LINE,O3P5_LINE,O4O5_LINE,O4O6_LINE,
O4O7_LINE,O4P1_LINE,O4P2_LINE,O4P3_LINE,
O4P4_LINE,O4P5_LINE,O5O6_LINE,O5O7_LINE,
O5P1_LINE,O5P2_LINE,O5P3_LINE,O5P4_LINE,
O5P5_LINE,O6O7_LINE,O6P4_LINE,O6P5_LINE,
O7P4_LINE,O7P5_LINE,P1P2_LINE,P1P3_LINE, 
P1P4_LINE, P1P5_LINE,P2P3_LINE, P2P4_LINE,
P2P5_LINE,P3P4_LINE, P3P5_LINE]


#siegbahn_dict  = dict((k.lower(), v) for k, v in zip(siegbahn_name,siegbahn_list))
#reverse_siegbahn_dict  = dict((k, v) for k, v in zip(siegbahn_list,siegbahn_name))
#shell_dict = dict((k.lower(), v) for k, v in zip(bindingE,shell_list))


def get_element_xrf_lines(Z,lines='all',incident_energy=100.0,lowE=0.0,highE=100.0,norm_to_one=False,include_escape=True,detector_element='Si'):
    """

    Generate a list of lines that a given element will generate for energy range listed

    

    
    """
    result={}
    xrf_list=[]
    escape_list=[]

    #
    #  This isn't very elegant...
    # 
    if lines.upper() == 'K':
        linelist = siegbahn_k_list
        namelist = siegbahn_k_name_list
    elif lines.upper() == 'KA':
        linelist = siegbahn_ka_list
        namelist = siegbahn_ka_name_list
    if lines.upper() == 'KB':
        linelist = siegbahn_kb_list
        namelist = siegbahn_kb_name_list
    elif lines.upper() == 'L':
        linelist =  siegbahn_l_list  
        namelist = siegbahn_l_name_list
    elif lines.upper() == 'LA' or lines.upper() == 'L1':
        linelist =  siegbahn_l1_list   
        namelist = siegbahn_l1_name_list
    elif lines.upper() == 'LB' or lines.upper() == 'L2':
        linelist =  siegbahn_l2_list   
        namelist = siegbahn_l2_name_list
    elif lines.upper() == 'LG' or lines.upper() == 'L3':
        linelist =  siegbahn_l3_list   
        namelist = siegbahn_l3_name_list
    elif lines.upper() == 'M':
        linelist =  siegbahn_m_list   
        namelist = siegbahn_m_name_list
    elif lines.upper() == 'MA':
        linelist =  siegbahn_ma_list   
        namelist = siegbahn_ma_name_list
    elif lines.upper() == 'MB':
        linelist =  siegbahn_mb_list   
        namelist = siegbahn_mb_name_list

    elif lines.upper() == 'MG':
        linelist =  siegbahn_mg_list   
        namelist = siegbahn_mg_name_list
    else:
        linelist = siegbahn_all_list     
        namelist = siegbahn_all_name_list
    
    
    for label,line in zip(namelist,linelist):
        lineE = xraylib.LineEnergy(Z,line)
        if lineE > 2.0 and lineE<=highE:    
            cs = xraylib.CS_FluorLine_Kissel(Z,line,incident_energy) 
            if cs > 0.0:
              #  xrf_list.append(XrayLine(label,lineE,cs))
                xrf_list.append([label,lineE,cs])
                if(include_escape):
                    ratio             = calc_escape_peak_ratios(lineE)
                    escape_energy     = calc_escape_peak_energy(lineE,detector_element)
                #    escape_list.append(EscapeLine(label,lineE,escape_energy,cs*ratio))
                    escape_list.append([label,lineE,escape_energy,cs*ratio])
 #   if norm_to_one:
  #     for line in xrf_list:
   #         line1[2] = line[2]/nsum
    #   for line in escape_list:
     #       line1[3] = line[3]/nsum
            
        
    result["lines"]=xrf_list
    result["escape"]=escape_list 
    return result

def totalxsecRange(element,energy_eV):
    
    
def 

def calc_filter_efficiency(paramdict,energy):
    """
    
    Based on the stuff between the sample and the detector
    work out how the XRF peaks should be attenuated...

    Three factors
    One is the sample attenuation of the x-rays...or matrix effect
    Next is the stuff between the sample and detector (air,detector windows etc.)
    Finally - is Enhancement factor. A major elemental component exciting a minor
    component which has a lower exciation energy.  

    mass attenuation, cross sections for a compound are calculated based on weighted
    averages.    cm2/g  - 
        

    Returns:  
    General attenuation curve - corresponding to air-detector losses
    Element specific factor corresponding to matrix effects
    
    The corrections are nicely outlined in...
    Fei He and Pierre J Van Espne
    ANALYTICAL CHEMISTRY, VOL. 63, NO, 20, OCTOBER 15, 1991 p 2237
    
    
    """
    
    #
    trans_curve=np.ones_like(energy)
    #
    # Step two  - the filters between the sample and detector
    # We assume these are normal to the beam...
    #
    pdict = paramdict["Acquisition_instrument"]["XRF"]
    if("Attenuators" in pdict):
        adict = pdict["Attenuators"]
        for attenuator in adict.keys():
                thickness = adict[attenuator]["thickness"]
                composition = adict[attenuator]["composition"]
                if(composition in paramdict["Materials"]):
                    density = paramdict["Materials"][composition]["Density"]
                    compoundList=paramdict["Materials"][composition]["CompoundList"]
                    fractionList=paramdict["Materials"][composition]["CompoundFraction"]
                else:
                    compoundList= composition
                    fractionList= 1.0
                    density = paramdict["Elements"][composition]["Density"]
                mass_thickness = density*thickness
                cross_section=calcTotalMassAttenuationCoefficients(compoundList,fractionList,energy,paramdict["Materials"],paramdict["Elements"])
                trans_curve = trans_curve * np.exp(-mass_thickness*cross_section)    
    return trans_curve


def calc_xrfdetector_efficiency(paramdict,energy):
    """
    
    Based on the stuff between the sample and the detector
    work out how the XRF peaks should be attenuated...

    Three factors
    One is the sample attenuation of the x-rays...or matrix effect
    Next is the stuff between the sample and detector (air,detector windows etc.)
    Finally - is Enhancement factor. A major elemental component exciting a minor
    component which has a lower exciation energy.  

    mass attenuation, cross sections for a compound are calculated based on weighted
    averages.    cm2/g  - 
        
    Returns:  
    General attenuation curve - corresponding to air-detector losses
    Element specific factor corresponding to matrix effects
    
    The corrections are nicely outlined in...
    Fei He and Pierre J Van Espne
    ANALYTICAL CHEMISTRY, VOL. 63, NO, 20, OCTOBER 15, 1991 p 2237
    
    
    """
    
    #
    trans_curve=np.ones_like(energy)
    detector = paramdict["Acquisition_instrument"]["XRF"]["Detector"]
    if("Attenuators" in detector):
        adict = detector["Attenuators"]
        for attenuator in adict.keys():
            thickness = adict[attenuator]["thickness"]
            composition = adict[attenuator]["composition"]
            if(composition in paramdict["Materials"]):
                density = paramdict["Materials"][composition]["Density"]
            else:
                density = paramdict["Elements"][composition]["Density"]
            # outgoing angle correction...
            # The thickness may be x microns but if x-rays aren't normal to it (on average) there will be some correction 
                
            mass_thickness = density*thickness
            cross_section=calcTotalMassAttenuationCoefficients(composition, 1.0 ,energy,paramdict["Materials"],paramdict["Elements"])
            if(attenuator!='sensor'):
                trans_curve = trans_curve * np.exp(-mass_thickness*cross_section)
            else:
                trans_curve = trans_curve * (1.0-np.exp(-mass_thickness*cross_section))
        
    return trans_curve



def calc_xrf_sample_matrix_correction(paramdict,energy):
    """
    
    Based on the stuff between the sample and the detector
    work out how the XRF peaks should be attenuated...

    Three factors
    One is the sample attenuation of the x-rays...or matrix effect
    Next is the stuff between the sample and detector (air,detector windows etc.)
    Finally - is Enhancement factor. A major elemental component exciting a minor
    component which has a lower exciation energy.  

    mass attenuation, cross sections for a compound are calculated based on weighted
    averages.    cm2/g  - 
        
    Returns:  
    General attenuation curve - corresponding to air-detector losses
    Element specific factor corresponding to matrix effects
    
    The corrections are nicely outlined in...
    Fei He and Pierre J Van Espne
    ANALYTICAL CHEMISTRY, VOL. 63, NO, 20, OCTOBER 15, 1991 p 2237
    
    
    """
    
    #
    #
    # Step three Matrix effects...
    #
    # The matrix correction is:
    # (1- exp(-chie()*density*thickness)/chie()
    #
    # where chie is:
    # cross_section(E_in)/incident_angle + cross_section(element line)/exit_angle
    # Loop over elements...and lines present...
    #
    absorption_correction = np.ones_like(energy)
    if("Matrix" in paramdict["Sample"].keys()):
        if paramdict["Sample"]["Matrix"] and "composition" in paramdict["Sample"]["Matrix"].keys():
            composition = paramdict["Sample"]["Matrix"]["composition"]
            if(composition in paramdict["Materials"].keys()):
                density = paramdict["Materials"][composition]["Density"]
                frac= paramdict["Materials"][composition]["CompoundFraction"]
                comp= paramdict["Materials"][composition]["CompoundList"]
            else:
                density = paramdict["Elements"][composition]["Density"]
                frac= 1.0
                comp= composition
            thickness = paramdict["Sample"]["Matrix"]["thickness"]
            incident_angle=math.cos(math.radians(paramdict["Acquisition_instrument"]["XRF"]["incident_angle"]))       
            exit_angle=math.cos(math.radians(paramdict["Acquisition_instrument"]["XRF"]["exit_angle"]))  
           
            mu_Ein  =calcTotalMassAttenuationCoefficients(comp,frac,paramdict["Acquisition_instrument"]["XRF"]["beam_energy"], paramdict["Materials"],paramdict["Elements"])
            mu_Eout =calcTotalMassAttenuationCoefficients(comp,frac,energy, paramdict["Materials"],paramdict["Elements"])
            mu_Tot = mu_Ein/incident_angle   + mu_Eout/exit_angle
            absorption_correction = (1.0-np.exp(-mu_Tot*density*thickness))/(mu_Tot*density*thickness)
        
    return absorption_correction

def calculate_transmission_corrections(paramdict,energy):
    trans = calc_filter_efficiency(paramdict,energy)
    trans *= calc_xrfdetector_efficiency(paramdict,energy)
    trans *= calc_xrf_sample_matrix_correction(paramdict,energy)
    return trans

def calcTotalMassAttenuationCoefficients(compoundList0, fractionList0,energy,Material,Element):
    """
    Usage:
        calcMassAttenuationCoefficients(compoundList, fractionList,
                                     energy = None,massfractions=False)
    compoundList - List of compounds into the material
    energy       - Energy at which the values are desired
    Material  : dictionary of pre-defined material properties
    Element   : dictionary of element properties
   
    The cross section is calculated using x-raylib.
    For a compound the total is the 
    weighted sum.
    u/p  =  sum_over_i (w_i * (u/p)i
    where w_i is the mass fraction of the ith component. 

    If the compound is a mixture of compounds
    e.g. compound = 0.4 SiO2 and 0.6 Fe
    then an average over those compounds will be calculated....

     
    """
    def __materialInCompoundList(lst,Material):
        for item in lst:
            if item in Material.keys():
                return True
        return False

    def calc_cross_section(z,energy):
        """
        Calculate the total cross section in cm^2/g
        using xraylib...
        Returns a numpy array...
        """
        result=np.zeros_like(energy)
        for i,en in np.ndenumerate(energy):
            result[i]=xraylib.CS_Total(z,en)
        return result


    if type(compoundList0) != type([]):
        compoundList = [compoundList0]
    else:
        compoundList = compoundList0
    if type(fractionList0) == np.ndarray:
        fractionList = fractionList0.tolist()
    elif type(fractionList0) != type([]):
        fractionList = [fractionList0]
    else:
        fractionList = fractionList0
    fractionList = [float(x) for x in fractionList]

    #
    # if one of the items listed in the compound 
    # is in the MATERIALS.DICT then we need to work out 
    # the average mass fraction for the total
    # E.g. if we have something with 0.4 SiO2 and 0.6 Si then we need
    # to work out the overall mass fraction of Si and O to calculate
    #  
    #

    while __materialInCompoundList(compoundList,Material):
        #
        # For each material in the compound...
        #
        for compound in compoundList:
            #
            # if it is a pre-defined materials
            #
            total=sum(fractionList)
            #allow materials in compoundList
            newcompound = []
            newfraction = []
            deleteitems = []
            
            if compound in Material.keys():
                if type(Material[compound]['CompoundList']) != type([]):
                    Material[compound]['CompoundList']=[Material[compound]['CompoundList']]
                if type(Material[compound]['CompoundFraction']) != type([]):
                    Material[compound]['CompoundFraction']=[Material[compound]['CompoundFraction']]
                Material[compound]['CompoundFraction'] = [float(x) for x in Material[compound]['CompoundFraction']]
                total = sum(Material[compound]['CompoundFraction'])
                j = compoundList.index(compound)
                compoundfraction = fractionList[j]
                i = 0
                for item in Material[compound]['CompoundList']:
                    newcompound.append(item)
                    newfraction.append(Material[compound]['CompoundFraction'][i] * compoundfraction /total)
                    i += 1
                deleteitems.append(j)
        if len(deleteitems):
            deleteitems.reverse()
            for i in deleteitems:
                del compoundList[i]
                del fractionList[i]
            for i in range(len(newcompound)):
                compoundList.append(newcompound[i])
                fractionList.append(newfraction[i])
    total=sum(fractionList)
    compoundFractionList = [float(x)/total for x in fractionList]
    materialElements = {}


    for compound in compoundList:
        elts=[]
        #get energy list
        if compound in Element.keys():
            elts=[compound]
            nbs =[1]
        else:
            elts= [ w for w in re.split('[0-9]', compound) if w != '' ]
            try:
                nbs= [ int(w) for w in re.split('[a-zA-Z]', compound) if w != '' ]
            except:
                raise ValueError("Compound '%s' not understood" % compound)
            if len(elts)==1 and len(nbs)==0:
                elts=[compound]
                nbs =[1]
        if (len(elts)==0 and len(nbs)==0) or (len(elts) != len(nbs)):
            print("compound %s not understood" % compound)
            raise ValueError("compound %s not understood" % compound)

        #the proportion of the element in that compound times the compound fraction
        fraction = [Element[elt]['Mass'] *nb for (elt, nb) in zip(elts, nbs) ]
        div      = sum(fraction)/compoundFractionList[compoundList.index(compound)]
        fraction = [x/div for x in fraction]

        for ele in elts:
            if ele not in materialElements.keys():
                materialElements[ele]  = fraction[elts.index(ele)]
            else:
                materialElements[ele] += fraction[elts.index(ele)]

    cs_total = np.ones_like(energy)

    for ele in materialElements.keys():
        z= Element[ele]['Z']
        cs_total=cs_total*calc_cross_section(z,energy)* materialElements[ele]
        
    return cs_total    

