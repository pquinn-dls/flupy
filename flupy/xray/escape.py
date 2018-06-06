import xraylib as xl
from math import log


def calc_escape_peak_energy(lineEnergy,detectortype):
    '''Escape peak energies...
     For Silicon you get a peak at E-1.742 (keV)
     For Ge there is a Ka and Kb at 9.876 and 10.984 (keV)
     so we get two offsets which need to be 
     weighted for Ka,Kb ratio 
     Parameters:
         
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


def calc_escape_peak_ratios(lineE,detectorelement='Si'):
    """
    Calculate ratio of escape peak to main peak based on emperical calculation
    from Alves et. al. 

    Parameters
    ----------
    detectorelement : string
        "Si" or "Ge"
 
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
        jump = xl.JumpFactor(14,xl.K_SHELL)
        fluy = xl.FluorYield(14,xl.K_SHELL)
        corr = fluy*(jump-1.0)/jump
        corr_photo = xl.CS_Photo(14,lineE)/xl.CS_Total(14,lineE)
        corr_trans = xl.RadRate(14,xl.KA_LINE)+xl.RadRate(14,xl.KB_LINE)
        mu_si= xl.CS_Total(14,lineE)
        mu_internal = xl.CS_Total(14,1.73998)
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
        if(lineE < 11.5):
            return 0.0,0.0
        jump = xl.JumpFactor(32,xl.K_SHELL)
        fluy = xl.FluorYield(32,xl.K_SHELL)
        corr = fluy*(jump-1.0)/jump
        corr_photo = xl.CS_Photo(32,lineE)/xl.CS_Total(32,lineE)
        corr_trans_ka = xl.RadRate(32,xl.KA_LINE)
        corr_trans_kb  =xl.RadRate(32,xl.KB_LINE)
        mu_ge= xl.CS_Total(32,lineE)#
        # one for the Ka and one for the Kb peak...
        mu_internal_ka = xl.CS_Total(32,xl.LineEnergy(32,xl.KA_LINE))
        r_ka = mu_internal_ka/mu_ge
        eta_ka = corr_trans_ka*corr_photo*corr*0.5*(1.0-r_ka*log(1.0+1.0/r_ka))
        ratio_ka =  eta_ka/(1.0-eta_ka)

        mu_internal_kb = xl.CS_Total(32,xl.LineEnergy(32,xl.KB_LINE))
        r_kb = mu_internal_kb/mu_ge
        eta_kb = corr_trans_kb*corr_photo*corr*0.5*(1.0-r_kb*log(1.0+1.0/r_kb))
        ratio_kb =  eta_kb/(1.0-eta_kb)

        return ratio_ka,ratio_kb
