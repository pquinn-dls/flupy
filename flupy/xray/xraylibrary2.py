# -*- coding: utf-8 -*-
"""
Created on Wed Apr 06 16:41:00 2011

@author: pq67
"""
import six
import xraylib
from flupy.xray.xrayliblines import xrayliblinelist,LinePair,\
                                          LineGroup,xraylibshelllist
from flupy.io.dict_io import loadFileFromDefault
from math import log
from collections import namedtuple,OrderedDict
import numpy as np
import re
import math
from itertools import product,combinations_with_replacement
import logging
from operator import add,mul
from functools import reduce

eV2keV = 1000.
sigma2fwhm = 2 * math.sqrt(2 * math.log(2))

xraylib.XRayInit()
xraylib.SetErrorMessages(0)
    
logger = logging.getLogger(__name__)

#
#  XRF model simulation and fitting is built around dictionarys containing 
#  these namedtupes  which contain the line name, location and intensity (cross section)
#
XrayLine = \
namedtuple('XrayLine', ('label','line_energy', 'intensity'))
EscapeLine = \
namedtuple('EscapeLine', ('label','line_energy','escape_energy', 'intensity'))
PileupLine = \
namedtuple('PileupLine', ('label','sublabels','line_energy', 'intensity'))
ComptonLine = \
namedtuple('ComptonLine', ('label','line_energy','intensity'))
ElasticLine = \
namedtuple('ElasticLine', ('label','line_energy','intensity'))


logger = logging.getLogger(__name__)

elementDB = loadFileFromDefault("ELEMENTS.DICT")
#print elementDB
elementDB = elementDB["Elements"]
# update the dictionary so you can reference by Z or full name also
elementDB.update({elm["Z"]: elm for elm in six.itervalues(elementDB)})
elementDB.update({elm["Name"].lower(): \
                  elm for elm in six.itervalues(elementDB)})


    

def cs_total(element,energy):
    """
    Get the total cross section for a given element when
    excitated with x-rays of a given energy (keV) 
    using the xraylib backend
    This is the total absorption cross section which is the 
    sum of the photoionization cross section, 
    the Rayleigh scattering cross section 
    and the Compton scattering cross section
    
    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar, list, tuple or numpy array
                    
    Returns:
        cs_total: float or numpy array - same size as input energy 
                  total cross section for this element at this incident 
                  energy in cm2/g.
    
    """
    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        xsec = np.zeros(len(energy))
        for i,enrg in enumerate(energy):
            xsec[i] = xraylib.CS_Total(z,enrg)
    else:
        xsec = xraylib.CS_Total(z,energy)
    return xsec

def cs_photo(element,energy):
    """
    Get the photo ionization cross section for a given element when
    excitated with x-rays of a given energy (keV) 
    using the xraylib backend
    
    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar, list, tuple or numpy array
                    
    Returns:
        cs_photo: float or numpy array - same size as input energy 
                  photoionization cross section for this element
                  at this incident energy in cm2/g.
    
    """
    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        xsec = np.zeroslike(energy)
        for i,enrg in enumerate(energy):
            xsec[i] = xraylib.CS_Photo(z,enrg)
    else:
        xsec = xraylib.CS_Photo(z,energy)
    return xsec

def cs_rayl(element,energy):
    """
    Get the rayleigh ionization cross section for a given element when
    excitated with x-rays of a given energy (keV) 
    using the xraylib backend
    
    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar, list, tuple or numpy array
                    
    Returns:
        cs_rayl: float or numpy array - same size as input energy 
                 rayleigh cross section for this element at this incident 
                 energy in cm2/g.
    
    """
    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        xsec = np.zeroslike(energy)
        for i,enrg in enumerate(energy):
            xsec[i] = xraylib.CS_Rayl(z,enrg)
    else:
        xsec = xraylib.CS_Rayl(z,energy)
    return xsec

def cs_compt(element,energy):
    """
    Get the compton ionization cross section for a given element when
    excitated with x-rays of a given energy (keV) 
    using the xraylib backend
    
    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar, list, tuple or numpy array
                    
    Returns:
        cs_compt: float or numpy array - same size as input energy 
                 rayleigh cross section for this element at this incident 
                 energy in cm2/g.
    
    """
    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        xsec = np.zeroslike(energy)
        for i,enrg in enumerate(energy):
            xsec[i] = xraylib.CS_Compt(z,enrg)
    else:
        xsec = xraylib.CS_Compt(z,energy)
    return xsec


def cs_total_kissel(element,energy):
    """
    """
    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        xsec = np.zeroslike(energy)
        for i,enrg in enumerate(energy):
            xsec[i] = xraylib.CS_Total_Kissel(z,enrg)
    else:
        xsec = xraylib.CS_Total_Kissel(z,energy)
    return xsec
    


def atomicWeight(element,Z):
    z = elementDB[element]["Z"]
    return xraylib.AtomicWeight(z)


def _lookupxlsubline(line):
    #
    #  This isn't very elegant...
    # 
    for xlline in xrayliblinelist:
        if xlline.label==line.lower():
            if isinstance(xlline,LinePair):
                return xlline
            elif isinstance(xlline,LineGroup):
                # if you've listed a group then just pick the first line
                print line," is a group of lines. Calculate using just",
                xlline.pairs[0].label
                return LinePair(xlline.pairs[0].sublabel,\
                                xlline.pairs[0].sublabel,\
                                xlline.pairs[0].subline)                        

def _lookupxlmainline(line):
    #
    #  This isn't very elegant...
    # 
    for xlline in xrayliblinelist:
        if xlline.label==line.lower():
            if isinstance(xlline,LinePair):
                return LineGroup(xlline.sublabel,(xlline))
            elif isinstance(xlline,LineGroup):
                return xlline

def _lookupxlshell(shell):
    
    for xlshell in xraylibshelllist:
        if xlshell.label==shell.upper():
           return xlshell

   
def _split_elementlinename(elementline):
    """
    Helper method 
    Splits the elementline input type  into
    element name and line
    E.g. Fe_K   is split into Fe, K
    Fe is split into Fe, All (all lines)
    
    Line names can be:
    
    Parameters:
        elementline :  String e.g.  Fe_Ka, Fe_Ka1 or Fe_K
        
    Returns:
        element name : element (string)
        line name    : line (string)

    """
    el = elementline.split('_')
    if len(el)!=2:
        element = el[0]
        lines = 'all'
    else:
        element = el[0]
        lines = el[1]
    return element,lines          

def line_energy(element,line):
    """
    
    Get the emission energy for a given element and line
    using the xraylib backend
    
    Parameters:
        Element : Name, Symbol or Atomic Number (Z) as input
        Line    : Accepts Siegbahn and line notations
                : line can be one of the following
                
                For line averages:
                'ka','lb','l1','l2','l3','la','lb','lg',
                
                or for the exact lines
                
                'ka1', 'ka2', 'ka3','kb1', 'kb2', 'kb3','kb4','kb5',
                'la1', 'la2', 'lb1', 'lb2', 'lb3', 'lb4', 'lb5',
                'lg1', 'lg2', 'lg3', 'lg4', 'll', 'ln',
                'ma1','ma2','mb','mg'
                    
    Returns:
        line energy: energy at which this fluorescence line appears 
    
    """
    z = elementDB[element]["Z"]
    if not isinstance(line,LinePair):
        line = _lookupxlsubline(line)
    return xraylib.LineEnergy(z,line.subline)                

def fluor_yield(element,shell):
    """
    
    Get the fluorescencec yield for a given element and shell
    using the xraylib backend
    
    Parameters:
        Element : Name, Symbol or Atomic Number (Z) as input
        Line    : Accepts Siegbahn and line notations
                : line can be one of the following
                
                For line averages:
                'ka','lb','l1','l2','l3','la','lb','lg',
                
                or for the exact lines
                
                'ka1', 'ka2', 'ka3','kb1', 'kb2', 'kb3','kb4','kb5',
                'la1', 'la2', 'lb1', 'lb2', 'lb3', 'lb4', 'lb5',
                'lg1', 'lg2', 'lg3', 'lg4', 'll', 'ln',
                'ma1','ma2','mb','mg'
                    
    Returns:
        line energy: energy at which this fluorescence line appears 
    
    """

    z = elementDB[element]["Z"]
    xlshell = _lookupxlshell(shell)
    return xraylib.FluorYield(z,xlshell.shell)        

def edge_energy(element,shell):
    """
    
    Get the absorption edge energy for a given element and shell
    using the xraylib backend
    
    Parameters:
        Element : Name, Symbol or Atomic Number (Z) as input
        Line    :'K','L1','L2','L3','La','Lb','Lg','Ma','Mb','Mg'
                    
    Returns:
        edge energy:   float 
    
    """

    z = elementDB[element]["Z"]
    xlshell = _lookupxlshell(shell)
    return xraylib.EdgeEnergy(z,xlshell.shell)    

def jump_factor(element, shell):
    """
    
    Get the absorption edge energy for a given element and shell
    using the xraylib backend
    
    Parameters:
        Element : Name, Symbol or Atomic Number (Z) as input
        Line    :'K','L1','L2','L3','La','Lb','Lg','Ma','Mb','Mg'
                    
    Returns:
        edge energy:   float 
    
    """

    z = elementDB[element]["Z"]
    xlshell = _lookupxlshell(shell)
    return xraylib.JumpFactor(z,xlshell.shell)

def cs_fluorline(element, line,excitation_energy):
    z = elementDB[element]["Z"]
    if not isinstance(line,LinePair):
        line = _lookupxlsubline(line)
    return xraylib.CS_FluorLine(z,line.subline,excitation_energy)
    
def cs_fluorline_kissel(element, line,excitation_energy):
    z = elementDB[element]["Z"]
    if not isinstance(line,LinePair):
        line = _lookupxlsubline(line)
    return xraylib.CS_FluorLine_Kissel(z,line.subline,excitation_energy)

def radrate(element,line):
    """
    
    Get the Fractional Radiative rate for a line using the xraylib backend
    
    Parameters:
        Element : Name, Symbol or Atomic Number (Z) as input
        Line    : Accepts Siegbahn and line notations
                : line can be one of the following
                
                For line averages:
                'ka','lb','l1','l2','l3','la','lb','lg',
                
                or for the exact lines
                
                'ka1', 'ka2', 'ka3','kb1', 'kb2', 'kb3','kb4','kb5',
                'la1', 'la2', 'lb1', 'lb2', 'lb3', 'lb4', 'lb5',
                'lg1', 'lg2', 'lg3', 'lg4', 'll', 'ln',
                'ma1','ma2','mb','mg'
                    
    Returns:
        radrate  :  the Fractional Radiative rate 
    
    """
    z = elementDB[element]["Z"]
    if not isinstance(line,LinePair):
        line = _lookupxlsubline(line)
    return xraylib.RadRate(z,line.subline)
  

def f1(element, energy):
    """
    
    The real part of anomalous X-ray scattering factor for 
    selected input energy (or energies) in eV.
    using xraylib backend

    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar, list, tuple or numpy array
                    
    Returns:
        f1 : float or numpy array - same size as input energy 
             real part of anomalous X-ray scattering factor at input
             energy or energies 
    
    """
    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        f1 = np.zeroslike(energy)
        for i,enrg in enumerate(energy):
            f1[i] = xraylib.Fi(z,enrg)
    else:
        f1 = xraylib.Fi(z,enrg)

    return f1

def f2(element, energy):
    """
    
    The imaginery part of anomalous X-ray scattering factor for 
    selected input energy (or energies) in eV.
    using xraylib backend

    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar, list, tuple or numpy array
                    
    Returns:
        f2 : float or numpy array - same size as input energy 
             imaginery part of anomalous X-ray scattering factor at input
             energy or energies 
    
    """

    z = elementDB[element]["Z"]
    if isinstance(energy, (list, tuple, np.ndarray)):
        f2 = np.zeroslike(energy)
        for i,enrg in enumerate(energy):
            f2[i] = xraylib.Fii(z,enrg)
    else:
        f2 = xraylib.Fii(z,enrg)
    return f2


def element_xrf_lines(element,incident_energy,lowE,highE,
                          exptdesc = None,
                          include_escape=False,
                          norm_to_one = False,
                          detector_element='Si'):

    """
    
    Generate a list of lines excited in element by this incident energy
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
    
    Parameters:
        element : Name, Symbol or Atomic Number (Z) as input
        energy  : energy in keV  - scalar
        lines   : One of the following options 'all' or a specific shell
        
                    
    Returns:
        lines : a list of named tuples 
                The named tuple has names     
                
    
    """

    xrf_list=[]
    # split the name if a line or shell is included e.g. Fe_K or Fe_Ka
    # if it's just an element symbol, e.g. Fe then all lines are considered...
    elementname,lines = _split_elementlinename(element)
    # Get all the xraylib macro lines for this line
    linegroups =  _lookupxlmainline(lines)
    # loop over the lines....
    for linep in linegroups.pairs:
        lineE   = line_energy(elementname,linep)
        if(lineE < lowE or lineE > highE):
            continue
        label   = linep.label
        cs = cs_fluorline_kissel(elementname,linep,incident_energy) 
        if exptdesc:
            cs = cs*expt_absorption_correction(exptdesc,
                                             lineE,incident_energy)
        xrf_list.append(XrayLine(elementname+"_"+label,lineE,cs))
        if cs > 0.0:
            if(include_escape):
                ratio             = calc_escape_peak_ratios(lineE)
                escape_energy     = calc_escape_peak_energy(lineE, \
                                                            detector_element)
                if escape_energy > 0.0:
                    xrf_list.append(EscapeLine(elementname+"_"+label+"_escape",\
                                               lineE,escape_energy,cs*ratio))
    if norm_to_one:
        sum_cs = sum([line.intensity \
                      for line in xrf_list if not isinstance(line,EscapeLine)])
        for i,line in enumerate(xrf_list):
            xrf_list[i] = line._replace(intensity=line.intensity/sum_cs)
    return xrf_list


def generate_element_lines(elementlist,incident_energy,
                          include_escape=False,
                          pileup_order = 2,
                          lowE = 0.0,
                          highE = np.inf,
                          detector_element='Si'):
    
    linedict = OrderedDict()
    for elementitem in elementlist:
        validlinelist =element_xrf_lines(elementitem,incident_energy,
                  lowE,highE,include_escape,
                  detector_element)
        if validlinelist:
            linedict[elementitem]=validlinelist

    return linedict

def generate_scatter_lines(incident_energy,scattering_angle,lowE,highE):
    linedict={}
    if(incident_energy< lowE or incident_energy > highE):
        return linedict 
    else:
        linedict["Elastic"]=[
                ElasticLine("elastic",incident_energy,1.0)]
        
        compton_e = compton_energy(incident_energy,scattering_angle)
        linedict["Compton"]=[ComptonLine("compton",compton_e,1.0)]
    return linedict

def generate_pileup_lines(linedict,pileup_order=2,verbose=False):
    """
    
    Generate the pileup lines from a starting dict containing one
    or more elements, compton or elastic lines. The dict will have been
    created by generate_detector_lines and/or generate_scatter_lines
    
    Pileup occurs when multiple photons hit the detector but the events
    are sp close together in time or at such a rate that the detector electronics
    sees them as single event. The result is that when the energy of this event 
    is recorded it will be the sum of the multiple photons.
    
    We can account for pileup peaks by just going through
    the possible permutations
    Example:
    if you have Fe_K and an Elastic line in the linedict
    then you will possibly have: 

    Fe_K+Fe_K
    Fe_K+Elastic
    Elastic+Elastic
    
    For a triple events the possible combinations would be..
    
    Fe_K+Fe_K+Fe_K
    Fe_K+Fe_K+Elastic
    Fe_K+Elastic+Elastic
    Elastic+Elastic+Elastic
    
    
    
    
    
    """
    new_linedict=OrderedDict()
    line_list =[key for key,value in linedict.iteritems() \
                if not isinstance(value,EscapeLine)]
    # loop from pair events up to triple, quadruple, or whatever 
    # pileup_order is set to.
    for nreps in range(2,pileup_order+1,1):
        # combinations of line names (e.g. Fe_k, Elastic)
        for pileup_combo in combinations_with_replacement(line_list,nreps):
            if verbose:
                print "***********************************************"
                print "pileup_combination",pileup_combo
                print "***********************************************"
            pileup_list=[]
            # From the combination of names extract the lines, excluding escape
            # peaks as these do not pileup
            llist=[[lline for lline in linedict[pc] \
                    if not isinstance(lline,EscapeLine)]\
                    for pc in pileup_combo]
            # for each line name you know have a list of lines
            # use product to get the combinations and then loop over that list
            for pp in product(*llist):
                if verbose:
                    print "pileup pairs",pp
                # sum the energy of the lines
                lineE =  reduce(add,[sline.line_energy for sline in pp])
                # multiply the intensity or probability of the lines
                intensity = reduce(mul,[sline.intensity for sline in pp])
                # create a list of the line labels and join them together
                sublabel  = [sline.label for sline in pp]
                label = '+'.join(sublabel)
                # Create a PileupLine namded tuple and add it to the list
                # of pileup lines for this particular combinaion 
                pileup_list.append(PileupLine(label,sublabel,
                                              lineE,intensity))
                if verbose:
                    print "Resulting pileup line",pileup_list[-1]
            # 
            new_linedict['+'.join(pileup_combo)] = pileup_list
    return new_linedict

#def filter_linedict(linedict,lowlim,highlim):
    
    


def calc_escape_peak_energy(lineEnergy,detectortype="Si"):
    '''
    Escape peak energies...
     For Silicon you get a peak at E-1.742 (keV)
     For Ge there is a Ka and Kb at 9.876 and 10.984 (keV)
     so we get two offsets which need to be 
     weighted for Ka,Kb ratio 
         
    Parameters
    ----------
         lineEnergy (keV)
         Detectortype  (string - 'Si', 'Si(Li)' or 'Ge') 
         if detector not recognisted it just returns the lineEnergy
     
    Returns
    -------
         List containing escape peak positions for selected detector
     
     '''
    
    if(detectortype=='Si' or detectortype=='Si(Li)'):
        return lineEnergy - 1.73998
    elif(detectortype == 'Ge'):
        return lineEnergy - 10.984,lineEnergy - 9.876
    else:
        return lineEnergy


def calc_escape_peak_ratios(lineEnergy,detectorelement='Si'):
    """
    Calculate ratio of escape peak to main peak based on emperical calculation
    from Alves et. al. 

    Parameters
    ----------
    lineEnergy      ; energy (keV)
    detectorelement : string "Si" or "Ge"
 
    Returns
    -------
    ratio of a single Si K escape peak to a single input peak 
    or
    ratio of Ge Ka, Kb escape peaks to a single input peak

    References
    ----------
    [1]
    "Experimental X-Ray peak-shape determination for a Si(Li) detector",
     L.C. Alves et al., Nucl. Instr. and Meth. in Phys. Res. B 109|110
    (1996) 129-133.  
    
    """
    
    if(detectorelement=='Si'):
        #
        # For Si the K peak is 95% of the transition
        # and the photoionization to total cross section is ~ 95% 
        # Si escape peak is typically only 0.2-1% (bigger at lower energies) 
        #
        jump = xraylib.JumpFactor(14,xraylib.K_SHELL)
        fluy = xraylib.FluorYield(14,xraylib.K_SHELL)
        corr = fluy*(jump-1.0)/jump
        corr_photo = \
        xraylib.CS_Photo(14,lineEnergy)/xraylib.CS_Total(14,lineEnergy)
        corr_trans = xraylib.RadRate(14,xraylib.KA_LINE)\
                    +xraylib.RadRate(14,xraylib.KB_LINE)
        mu_si= xraylib.CS_Total(14,lineEnergy)
        mu_internal = xraylib.CS_Total(14,1.73998)
        r = mu_internal/mu_si
        eta = corr_trans*corr_photo*corr*0.5*(1.0-r*log(1.0+1.0/r))
        ratio =  eta/(1.0-eta)
        #
        # escape peak sigma should be narrower than  the main peak.
        #
        return ratio
    else:
        # 
        # Ge detector...
        # Ge has a large escape peak ratio ~ 5-15% and a Ka and kb component
        #
        if(lineEnergy < 11.5):
            return 0.0,0.0
        jump = xraylib.JumpFactor(32,xraylib.K_SHELL)
        fluy = xraylib.FluorYield(32,xraylib.K_SHELL)
        corr = fluy*(jump-1.0)/jump
        corr_photo = \
        xraylib.CS_Photo(32,lineEnergy)/xraylib.CS_Total(32,lineEnergy)
        corr_trans_ka = xraylib.RadRate(32,xraylib.KA_LINE)
        corr_trans_kb  =xraylib.RadRate(32,xraylib.KB_LINE)
        mu_ge= xraylib.CS_Total(32,lineEnergy)#
        # one for the Ka and one for the Kb peak...
        mu_internal_ka = \
            xraylib.CS_Total(32,xraylib.LineEnergy(32,xraylib.KA_LINE))
        r_ka = mu_internal_ka/mu_ge
        eta_ka = corr_trans_ka*corr_photo*corr*0.5*(1.0-r_ka*log(1.0+1.0/r_ka))
        ratio_ka =  eta_ka/(1.0-eta_ka)

        mu_internal_kb = \
        xraylib.CS_Total(32,xraylib.LineEnergy(32,xraylib.KB_LINE))
        r_kb = mu_internal_kb/mu_ge
        eta_kb = corr_trans_kb*corr_photo*corr*0.5*(1.0-r_kb*log(1.0+1.0/r_kb))
        ratio_kb =  eta_kb/(1.0-eta_kb)

        return ratio_ka,ratio_kb

def compton_energy(energy,scattering_angle):
    """
    
    The scattering of photons from charged particles 
    is called Compton scattering after Arthur Compton who was the first to 
    measure photon-electron scattering in 1922. When the incoming photon 
    gives part of its energy to the electron, then the scattered photon has 
    lower energy and according to the Planck relationship has lower 
    frequency and longer wavelength. The wavelength change in such 
    scattering depends only upon the angle of scattering for a 
    given target particle
    
    Parameters:
    Inputs:
        
        energy (keV)  - input photon energy
        scattering_angle  (degrees) - compton scattering angle 
    
    Output:
        
        compton energy (keV) - location of the compton peak        
    
    
    """

    mc2 = 511
    comp_denom = (1. + energy / mc2 *
                 (1. - np.cos(np.deg2rad(scattering_angle))))
    compton_peak = energy / comp_denom
	
    return compton_peak    



def calc_filter_efficiency(exptdesc,line_energy):
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
    # Step two  - the filters between the sample and detector
    # We assume these are normal to the beam...
    #
    if isinstance(line_energy, (list, tuple, np.ndarray)):
        trans_curve=np.ones_like(line_energy)
    else:
        trans_curve=1.0
    if isinstance(exptdesc, dict):        
        try:
            pdict   = exptdesc["Acquisition_instrument"]["XRF"]
            matdict = exptdesc["Materials"]
            eldict  = exptdesc["Elements"]
            if("Attenuators" in pdict):
                attdict = pdict["Attenuators"]
                for attenuator in attdict.keys():
                        thickness = attdict[attenuator]["thickness"]
                        composition = attdict[attenuator]["composition"]
                        if(composition in matdict):
                            density = matdict[composition]["Density"]
                            compoundList=matdict[composition]["CompoundList"]
                            fractionList=matdict[composition]["CompoundFraction"]
                        else:
                            compoundList= composition
                            fractionList= 1.0
                            density = eldict[composition]["Density"]
                        mass_thickness = density*thickness
                        cross_section=calcTotalMassAttenuationCoefficients(
                                compoundList,fractionList,line_energy,matdict,eldict)
                        trans_curve = trans_curve*np.exp(-mass_thickness*cross_section)    
            return trans_curve
        except KeyError:
            print "expt description not provided or incomplete \
            so no correction applied"
            return trans_curve
    else:
        return trans_curve


def calc_xrfdetector_efficiency(exptdesc,line_energy):
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
    
    if isinstance(line_energy, (list, tuple, np.ndarray)):
        trans_curve=np.ones_like(line_energy)
    else:
        trans_curve=1.0
    if isinstance(exptdesc, dict):
        try:
            detector = exptdesc["Acquisition_instrument"]["XRF"]["Detector"]
            if("Attenuators" in detector):
                adict = detector["Attenuators"]
                for attenuator in adict.keys():
                    thickness = adict[attenuator]["thickness"]
                    composition = adict[attenuator]["composition"]
                    if(composition in exptdesc["Materials"]):
                        density = exptdesc["Materials"][composition]["Density"]
                    else:
                        density = exptdesc["Elements"][composition]["Density"]
                    # outgoing angle correction...
                    # The thickness may be x microns but if x-rays aren't 
                    # normal to it (on average) there will be some correction 
                        
                    mass_thickness = density*thickness
                    cross_section=calcTotalMassAttenuationCoefficients(composition, 
                               1.0 ,line_energy,exptdesc["Materials"],
                               exptdesc["Elements"])
                    if(attenuator!='sensor'):
                        trans_curve = trans_curve * \
                            np.exp(-mass_thickness*cross_section)
                    else:
                        trans_curve = trans_curve * \
                            (1.0-np.exp(-mass_thickness*cross_section))
            return trans_curve
    
        except KeyError:
            print "Main expt keys not present...Materials,Elements,\
                Acquisition_instrument/XRF...returning ones "
            return trans_curve
            
    else:
        return trans_curve



def calc_xrf_sample_matrix_correction(exptdesc,line_energy,incident_energy):
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
    if isinstance(line_energy, (list, tuple, np.ndarray)):
        absorption_correction=np.ones_like(line_energy)
    else:
        absorption_correction=1.0
    if isinstance(exptdesc, dict):
        try:
            if("Matrix" in exptdesc["Sample"].keys()):
                if exptdesc["Sample"]["Matrix"] and \
                    "composition" in exptdesc["Sample"]["Matrix"].keys():
                    composition = exptdesc["Sample"]["Matrix"]["composition"]
                    # is it a defined material or element...
                    # extract the density and composition
                    if(composition in exptdesc["Materials"].keys()):
                        density = exptdesc["Materials"][composition]["Density"]
                        frac= exptdesc["Materials"][composition]["CompoundFraction"]
                        comp= exptdesc["Materials"][composition]["CompoundList"]
                    else:
                        density = exptdesc["Elements"][composition]["Density"]
                        frac= 1.0
                        comp= composition
                    # sample thickness
                    thickness = exptdesc["Sample"]["Matrix"]["thickness"]
                    # beam angles for path length calculations
                    incident_angle=math.cos(math.radians(\
                        exptdesc["Acquisition_instrument"]["XRF"]["incident_angle"]))       
                    exit_angle=math.cos(math.radians(\
                            exptdesc["Acquisition_instrument"]["XRF"]["exit_angle"]))  
                    # 
                    mu_Ein  =calcTotalMassAttenuationCoefficients(comp,frac,\
                        incident_energy, \
                        exptdesc["Materials"],exptdesc["Elements"])
                    mu_Eout =calcTotalMassAttenuationCoefficients(comp,frac,
                                                                  line_energy, 
                        exptdesc["Materials"],exptdesc["Elements"])
                    mu_Tot = mu_Ein/incident_angle   + mu_Eout/exit_angle
                    absorption_correction = (1.0-\
                        np.exp(-mu_Tot*density*thickness))/(mu_Tot*density*thickness)
            return absorption_correction
    
        except (TypeError,KeyError):
            print "Main expt keys not present...Materials,Elements,\
                Acquisition_instrument/XRF...returning ones "
            return absorption_correction
    else:
        return absorption_correction

def expt_absorption_correction(exptdesc,line_energy,incident_energy):
    """
    """
    # attenuation from filters between sample and detector
    trans = calc_filter_efficiency(exptdesc,line_energy)
    # attenuation/absorption efficiency of detector
    trans *= calc_xrfdetector_efficiency(exptdesc,line_energy)
    # sample attenuation - path length self-absorption
    trans *= calc_xrf_sample_matrix_correction(exptdesc,line_energy,incident_energy)
    return trans

def calcTotalMassAttenuationCoefficients(compoundList0, fractionList0,
                                         energy,Material,Element):
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
                    Material[compound]['CompoundList']=\
                        [Material[compound]['CompoundList']]
                if type(Material[compound]['CompoundFraction']) != type([]):
                    Material[compound]['CompoundFraction']=\
                        [Material[compound]['CompoundFraction']]
                Material[compound]['CompoundFraction'] = \
                    [float(x) for x in Material[compound]['CompoundFraction']]
                total = sum(Material[compound]['CompoundFraction'])
                j = compoundList.index(compound)
                compoundfraction = fractionList[j]
                i = 0
                for item in Material[compound]['CompoundList']:
                    newcompound.append(item)
                    newfraction.append(
                            Material[compound]['CompoundFraction'][i] * 
                            compoundfraction /total)
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
        cs_total=cs_total*cs_total(z,energy)* materialElements[ele]
        
    return cs_total    







def _get_element_and_line(xray_line):
    """
    Returns the element name and line character for a particular X-ray line as
    a tuple.

    By example, if xray_line = 'Mn_Ka' this function returns ('Mn', 'Ka')
    """
    lim = xray_line.find('_')
    if lim == -1:
        raise ValueError("Invalid xray-line: %" % xray_line)
    return xray_line[:lim], xray_line[lim + 1:]


def _get_energy_xray_line(xray_line):
    """
    Returns the energy (in keV) associated with a given X-ray line.

    By example, if xray_line = 'Mn_Ka' this function returns 5.8987
    """
    element, line = _get_element_and_line(xray_line)
    return elements_db[element]['Atomic_properties']['Xray_lines'][
        line]['energy (keV)']


def _get_xray_lines_family(xray_line):
    """
    Returns the family to which a particular X-ray line belongs.

    By example, if xray_line = 'Mn_Ka' this function returns 'Mn_K'
    """
    return xray_line[:xray_line.find('_') + 2]


def _parse_only_lines(only_lines):
    if isinstance(only_lines, str):
        pass
    elif hasattr(only_lines, '__iter__'):
        if any(isinstance(line, str) is False for line in only_lines):
            return only_lines
    else:
        return only_lines
    only_lines = list(only_lines)
    for only_line in only_lines:
        if only_line == 'a':
            only_lines.extend(['Ka', 'La', 'Ma'])
        elif only_line == 'b':
            only_lines.extend(['Kb', 'Lb1', 'Mb'])
    return only_lines


def get_xray_lines_near_energy(energy, width=0.2, only_lines=None):
    """Find xray lines near a specific energy, more specifically all xray lines
    that satisfy only_lines and are within the given energy window width around
    the passed energy.

    Parameters
    ----------
    energy : float
        Energy to search near in keV
    width : float
        Window width in keV around energy in which to find nearby energies,
        i.e. a value of 0.2 keV (the default) means to search +/- 0.1 keV.
    only_lines :
        If not None, only the given lines will be added (eg. ('a','Kb')).

    Returns
    -------
    List of xray-lines sorted by energy difference to given energy.
    """
    only_lines = _parse_only_lines(only_lines)
    valid_lines = []
    E_min, E_max = energy - width / 2., energy + width / 2.
    for element, el_props in elements_db.items():
        # Not all elements in the DB have the keys, so catch KeyErrors
        try:
            lines = el_props['Atomic_properties']['Xray_lines']
        except KeyError:
            continue
        for line, l_props in lines.items():
            if only_lines and line not in only_lines:
                continue
            line_energy = l_props['energy (keV)']
            if E_min <= line_energy <= E_max:
                # Store line in Element_Line format, and energy difference
                valid_lines.append((element + "_" + line,
                                    np.abs(line_energy - energy)))
    # Sort by energy difference, but return only the line names
    return [line for line, _ in sorted(valid_lines, key=lambda x: x[1])]


def get_FWHM_at_Energy(energy_resolution_MnKa, E):
    """Calculates an approximate FWHM, accounting for peak broadening due to the
    detector, for a peak at energy E given a known width at a reference energy.

    The factor 2.5 is a constant derived by Fiori & Newbury as references
    below.

    Parameters
    ----------
    energy_resolution_MnKa : float
        Energy resolution of Mn Ka in eV
    E : float
        Energy of the peak in keV

    Returns
    -------
    float : FWHM of the peak in keV

    Notes
    -----
    This method implements the equation derived by Fiori and Newbury as is
    documented in the following:

        Fiori, C. E., and Newbury, D. E. (1978). In SEM/1978/I, SEM, Inc.,
        AMF O'Hare, Illinois, p. 401.

        Goldstein et al. (2003). "Scanning Electron Microscopy & X-ray
        Microanalysis", Plenum, third edition, p 315.

    """
    FWHM_ref = energy_resolution_MnKa
    E_ref = _get_energy_xray_line('Mn_Ka')

    FWHM_e = 2.5 * (E - E_ref) * eV2keV + FWHM_ref * FWHM_ref

    return math.sqrt(FWHM_e) / 1000.  # In mrad

