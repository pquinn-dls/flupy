#
# imports
#
import _xraylib as xl
import numpy as np
from itertools import combinations_with_replacement,combinations

#
# Local imports
#
from flupy.algorithms.xrf_calculations.lineshapes import gaussian,gaussian_tail,shelf


def pileup_characteristic_lines(x, Z,translist,transitions, sigma, tail, slope, step,correction=None):
    """
    Calculates characteristic XRF lines for pileup...
    E.g. For K lines we need Ka+Ka, Ka+Kb and Kb+Kb
    the amplitudes will be the product of the transition rates...

    Args:
    -----

            x:              Array of energy in KeV.
            Z:              Input element number
            transitions:    List of macros for denoting transition
            sigma:          Detector precision.
            tail_amplitude: Amplitude for gaussian_tail contribution.
            slope:          Slope of tail. Slope = sigma.
            S:              Step.
            pileup          include pileup peak

    Returns:
    --------

            Returns characteristic XRF lines for given transitions in element Z.
    """

    tline  = np.zeros(x.size)
      
    list_of_trans=[]
    list_of_amps=[]
    for i,gt in enumerate(transitions):
        if(translist[i]!=0):
            for ht in gt:
                amplitude = xl.RadRate(Z, ht)
                list_of_amps.append(amplitude)
                list_of_trans.append(ht)
    index_c = range(len(list_of_trans))
    #
    # use itertools to work out combinations of transitions
    # e.g. (1,1),(1,2),(2,2) etc.
    # A bit lengthy but I think it get the job done correctly..
    #
    combinations=combinations_with_replacement(index_c,2)
    for comb in combinations:
        amp1 =list_of_amps[comb[0]] 
        amp2 =list_of_amps[comb[1]]
        if (amp1 != 0.0 and amp2 !=0.0):
            energy1 = xl.LineEnergy(Z,list_of_trans[comb[0]])
            energy2 = xl.LineEnergy(Z,list_of_trans[comb[1]])
            # amplitude is the ratio of the peak probabilities..
            amplitude = amp1*amp2
            # energy is the sum of these transition energys
            energy = energy1+energy2
            tline = tline + amplitude * (gaussian(x, energy, sigma)[0] +
                    gaussian_tail(x, energy, energy, sigma,tail, slope) +
                                 shelf(x, energy, sigma, step)) # //+

    if(correction!=None):
        tline=tline*correction        
    area=np.sum(tline)                                 

    line = tline
       
    return line,area



def sumup_line_pair(x, Z1, Z2, translist1, translist2,transitions, sigma, tail, slope, step,correction=None):
    """
    Calculates characteristic XRF lines for pileup...
    E.g. For K lines we need Ka+Ka, Ka+Kb and Kb+Kb
    the amplitudes will be the product of the transition rates...

    Args:
    -----

            x:              Array of energy in KeV.
            Z:              Input element number
            transitions:    List of macros for denoting transition
            sigma:          Detector precision.
            tail_amplitude: Amplitude for gaussian_tail contribution.
            slope:          Slope of tail. Slope = sigma.
            S:              Step.
            pileup          include pileup peak

    Returns:
    --------

            Returns characteristic XRF lines for given transitions in element Z.
    """

    tline  = np.zeros(x.size)
    validtrans1=[transitions[i] for i, v in enumerate(translist1) if v > 0]
    validtrans2=[transitions[i] for i, v in enumerate(translist2) if v > 0]
    validtrans =validtrans1+validtrans2 
    transind = range(len(validtrans))
    pileupcombos = combinations_with_replacement(transind, 2)
     
    for combo in pileupcombos:
        if combo[0]<len(validtrans1):
            za=Z1
            zb=Z2
        else:
            za=Z2
            zb=Z1
            
        energy1=0.0
        amplitude1=0.0
        energy2=0.0
        amplitude2=0.0
        validtrans[combo[0]] + validtrans[combo[1]]
        for ht in validtrans[combo[0]]:
            # rad rate for relative size of the transition
            amplitude = xl.RadRate(za, ht)
            if amplitude != 0.0:
                # energy of this transition..
                energy1 += amplitude*xl.LineEnergy(za, ht)
                amplitude1+=amplitude
        energy1=energy1/amplitude1
        for ht in validtrans[combo[1]]:
            # rad rate for relative size of the transition
            amplitude = xl.RadRate(zb, ht)
            if amplitude != 0.0:
                # energy of this transition..
                energy2 += amplitude*xl.LineEnergy(zb, ht)
                amplitude2+=amplitude
        energy2=energy2/amplitude2
        amplitude = amplitude1*amplitude2
        energy = energy1+energy2
         
        amplitude = xl.RadRate(Z1, ht)
        if combo[0]>=len(validtrans1):
                amplitude = xl.RadRate(Z2, ht)
        tline = tline + amplitude * (gaussian(x, energy, sigma)[0] +
                    gaussian_tail(x, energy, energy, sigma,tail, slope) +
                                 shelf(x, energy, sigma, step)) # //+
            
    if(correction!=None):
        tline=tline*correction        
    area=np.sum(tline)                                 
    line = tline
       
    return line,area