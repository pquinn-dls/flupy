
import numpy as np
from scipy.optimize import nnls


class linlsqconstrained2:
    
    def __init__(self, a,constraints,h=None):
        self.a=a
        self.constraints = constraints
        #self.u,self.n2,self.G_bar,self.vS_inv,self.fprime = \
        #self.precalculate_useful_matrices(a,constraints)
        
        if(h==None):
            self.h = np.zeros(constraints.shape[0])

    def precalculate_useful_matrices(self,DB, constraints,bp):
        # Step 1 - convert to an LDP problem
        mod_weight = 1. / np.sqrt(bp)
        rweight = np.diag(mod_weight)
        a  = np.dot(rweight.T, DB)
        
        u,s,vh = np.linalg.svd(a,full_matrices=True)
        v=vh.T
        n2 = a.shape[1]
        # Pre- calculate some of the arrays....
        S = np.zeros((a.shape[1],a.shape[1]))
        S[:,:] = np.diag(s)
        S_inv = np.linalg.inv(S)
        vS_inv = np.dot(v,S_inv)
        # big G    
        fprime = np.zeros(n2+1)
        fprime[-1] = 1.0
        
        G_bar = np.dot(constraints,vS_inv)
        
        return u,n2,G_bar,vS_inv,fprime 

    def solve_for_x(self,bp): 
        """
         solve for a given input spectrum...    
        """       
        #u,n2,h,G_bar,bp,fprime,vS_inv
        self.u,self.n2,self.G_bar,self.vS_inv,self.fprime = \
            self.precalculate_useful_matrices(self.a,self.constraints,bp)
        
        fp     = np.dot(self.u.T,bp)[:self.n2]
        offset = np.dot(self.G_bar,fp)
        h_bar  = self.h - offset
        ht     = h_bar.reshape((h_bar.shape[0],1))
        E      = np.vstack((self.G_bar.T,ht.T))
        mpp    = nnls(E,self.fprime)
        ep     = np.dot(E,mpp[0]) - self.fprime
        mp     = -ep[:self.n2]/ep[self.n2]
        mest   = np.dot(self.vS_inv ,fp+mp)
        return mest
# def linearlsq_with_inequality_constraints(a,b, constraints,h=None):
#     """
#     Solve Ax=B  subject to constraints*x>=h     
#     h is optional and if not included is set to zeros
#     
#     """
#     # Step 1 - convert to an LDP problem
#     if(h==None):
#         h = np.zeros(constraints.shape[0])
#     u,n2,G_bar,vS_inv,fprime = __precalculate_useful_matrices(a,constraints,h)
#     mest = __solve_for_x(u,n2,h,G_bar,b,fprime,vS_inv)
#     return mest
# 
#         
# def __solve_for_x(u,n2,h,G_bar,bp,fprime,vS_inv): 
#     """
#      solve for a given input spectrum...    
#     """       
#     fp     = np.dot(u.T,bp)[:n2]
#     offset = np.dot(G_bar,fp)
#     h_bar  = h - offset
#     ht     = h_bar.reshape((h_bar.shape[0],1))
#     E      = np.vstack((G_bar.T,ht.T))
#     mpp    = nnls(E,fprime)
#     ep     = np.dot(E,mpp[0]) - fprime
#     mp     = -ep[:n2]/ep[n2]
#     mest   = np.dot(vS_inv ,fp+mp)
#     return mest
#     
# 

