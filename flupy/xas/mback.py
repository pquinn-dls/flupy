
import math
import numpy
import scipy.optimize
import scipy.special
import xraylib
from xray import elam
import numpy as np

def mbgfun(alpha,energy,Mu_obs,E_edge,E_fluo,weight):
    
    scale= alpha[0]	# scale factor for Mu_obs
    amp  = alpha[1]	# Amplitude of erfc()
    xi   = alpha[2]	    # spectra width in erfc()
    coef = alpha[3:]	# coefficients for polynomials
	 #
	 # back ground - pre-edge
    #
    bg_poly=numpy.polyval(coef,energy-E_edge)
	 #
	 # erfc + background
	 #
    bg_erfc= amp*scipy.special.erfc((energy-E_fluo)/xi)
	 bg_tot = bg_erfc+bg_poly
	 #
	 #
	 #
    Mu=Scale*(Mu_obs-bg_tot)*weight    
    return Mu,bg_tot,bg_erfc,bg_poly
	

    
def mback_residual(alpha,energy,mu_tab,mu_obs,E_edge,E_fluo,weight):
    mu,back,back1,back2=mbgfun(alpha,energy,mu_obs,E_edge,E_fluo,weight)
    difference=mu-mu_tab
    #print numpy.sum(mu*mu),len(mu),alpha
    return difference

    
def erfc_residual(guess,energy,mu,E_edge,E_fluo):
    bg_erfc=guess[0]*scipy.special.erfc((energy-E_fluo)/guess[1])
    difference=bg_erfc-mu
    return difference

mbackxraydb = xraylibDB()

def mback(edge_energy,fluo_energy,energy,mu_raw,exclusion_low=20,exclusion_high=40,order=2):
    """
	
    MBACK Normalize raw data to tabulated mass absorption coefficients
    with erfc() and polynomials as the background functions.
    
    Parameters
    ----------
    element : string  E.G. "Fe" 
    edge    : string "K", "L3","L2","L1" etc. 
    energy  : array of energy values
    mu_raw  : the raw data intensity values
    fitrange: 
    order   : order of polynomial used in fit - default = 2    
 
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

    
    Converted to python - P Quinn 2010
	
    """    
    
    db=elam.elam()
    # (1a)  energy in KeV
    E_edge=db.edgeEnergy(element,edge)
    fluo=db.lines(element,edge)
    # Lets just pick the line with max intensity...for this edge)
    E_fluo=fluo[0][2]
    #for data in fluo:
    #    E_fluo=max(E_fluo,data[3])
    #
    print 'E_fluo',E_fluo,E_edge,energy_eV[0]
    #
    # pre-edge region, post-edge region and total fitted region
    #    
    pre_edge  = E_edge - abs(exclusion_low)
    post_edge = E_edge + abs(exclusion_high)
    idx_pre   = energy_eV<=pre_edge
    idx_post  = energy_eV>=post_edge
    idx_xanes = np.logical_or(energy_eV>pre_edge,energy_eV<post_edge)
    # Now extract the data defined by the fit-range - ignoring the XANES		
    idx_fit=np.logical_or(idx_pre,idx_post) # index for pre-/post-edge data points
    num_pre_edge_points = np.sum(idx_pre)
    num_post_edge_points= np.sum(idx_post)
    #
    # (2) Tabulated mass absorption coefficients
	 #
    mutab=db.totalxsecRange(element,energy_eV)
	 #
	
    #
    # Points below the edge.....
    # Fit the pre-edge and remove it
    #
    idx_res= energy_eV<E_edge
    # Fit 2nd order poly to region < edge and then get poly curve over the full energy range
    poly=numpy.polyfit(energy_eV[idx_res],mutab[idx_res],Nth)
    # interpolated over full energy range..
    res=numpy.polyval(poly,energy_eV)
    # mu_tabulated = mu_tabluted - background poly
    mutab = mutab-res
	 # 
    # Should be zero as we have removed the background but just set to zero....
    # (3a) Scale raw data and mu_tab
    # Max over the xanes range
    mutab[idx_res]=0.0
    mu_span = mutab.max()
    # upper and lower limits of the data over the xanes region
    mu_xanes = mu_raw[idx_xanes]
    raw_span=mu_xanes.max()-mu_xanes.min()
    # scale mu 
    # extract subset of fit data and divide by span
    mu_tab=mutab[idx_fit]/mu_span
    #mu_tab=mutab[idx_fit]
    # scale the raw data
    #mu_obs=mu_raw[idx_fit]/raw_span
    mu_obs=mu_raw[idx_fit]
    # energy fit data set
    eV_obs=energy_eV[idx_fit]
    # (3b) Normalization
    xi0=500
    # Fit a poly of order n from pre-edge to end  and use as your initial guess
    pi0=numpy.polyfit(eV_obs[num_pre_edge_points+1:]-E_edge,mu_obs[num_pre_edge_points+1:]-(mu_obs.max()-mu_obs.min()),Nth)
    #pi0=numpy.polyfit(KeV_obs[num_pre_edge_points+1:]-E_edge,mu_tab[num_pre_edge_points+1:],Nth)
    # Guess at the erfc function...
    #fits=scipy.optimize.leastsq(erfc_residual,[1,500.0],args=(eV_obs,mu_obs-mu_xanes.min(),E_edge,E_fluo),full_output=1)
    
    #print 'fits!!!',fits[0]
    # S,A,xi, pi0
    guess=numpy.array([1.0,10.0,xi0])
    guess=numpy.append(guess,pi0) # initial guess: [S,A,xi,{P_i}]
    # A list of ones
    weight_one=numpy.ones(len(mu_raw))
    # Weighting for the fit.............weight to number of points pre and post edge
    weight_fit=numpy.append(numpy.tile(1.0/math.sqrt(num_pre_edge_points),num_pre_edge_points),\
    numpy.tile(1.0/math.sqrt(num_post_edge_points),num_post_edge_points))
    
    fits=scipy.optimize.leastsq(mback_residual,guess,args=(eV_obs,mu_tab*weight_fit,mu_obs,E_edge,E_fluo,weight_fit),full_output=1)
    # Normalized data is with no weighting....
    print 'fit results',fits[0]
    #norm_xas,back,back1,back2=mbgfun(fits[0],energy_eV,mu_raw/raw_span,E_edge,E_fluo,weight_one)
    #return norm_xas,back*raw_span/fits[0][0],back1*raw_span/fits[0][0],back2*raw_span/fits[0][0],mutab
    norm_xas,back,back1,back2=mbgfun(fits[0],energy_eV,mu_raw,E_edge,E_fluo,weight_one)    
    return norm_xas,back,mutab/mu_span

	
#fin=open(r"c:\code\dataprocess\xas\1103.dat.xmu")
#data=fin.readlines()
#fin.close()
#enrg=[]
#mux=[]
#for i in range(1,len(data)):
#	spx=data[i].split()
	#print spx
	#enrgt=float(spx[0])
	#enrg.append(enrgt)
	#muxt=float(spx[25])/float(spx[3])
	#muxt=float(spx[1])
	#mux.append(muxt)
#result,back,mutab=mback("Ti","K",numpy.array(enrg),numpy.array(mux),fitrange=[15,50],Nth=2)
#fout=open(r"c:\code\dataprocess\xas\1103norm.dat","w")
#for i in range(len(result)):
#    print >> fout,enrg[i],result[i],mux[i],back[i],mutab[i]
#fout.close()    
    



