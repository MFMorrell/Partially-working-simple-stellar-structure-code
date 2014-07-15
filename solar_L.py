#This file contains functions for calculation of DL/DM and its derivatives. The reactions use to calculate energy  are the Proton Proton, CNO and Triple alpha .

#Import
from solar_rho import rho, rho_deriv
from MESA import Tsun, Msun, Lsun
from scipy.interpolate import splrep, splev

from numpy import exp
#constants
eppCon= 2.41e6 #ergs
ecnoCon=8.67e27 #ergs
e3aCon=5.09e11 #ergs

#Scaling factor
s = Msun/Lsun     
#
    

#Scaled Proton Proton reaction
def enucpp(P, T, M, X):
  
    #find rho
    rhof = rho(P, T, M, X)
    # find and interpolate mass fractions
    [fH, fHe, fOther] = X
    H = splev(M, fH)
    #transform T    
    T6=(T*Tsun)*1.0e-6
    #refactor and evaluate
    c = eppCon*H**2.0
    Texp = T6**(-2.0/3.0)*exp(-33.80*T6**(-1.0/3.0))
    
    return c*rhof*Texp*s
  
#Scaled Proton Proton reaction derivative  
def enucpp_deriv(P, T, M, X):
    
    #find rho
    rhof = rho(P, T, M, X)
    rhod = rho_deriv(P, T, M, X)
    # find and interpolate mass fractions
    [fH, fHe, fOther] = X
    H = splev(M, fH)
    #transform T 
    T6=(T*Tsun)*1.0e-6
    Tsun6 = Tsun*1.0e-6
    #refactor and evaluate
    c = eppCon*H**2.0
    Texp = T6**(-2.0/3.0)*exp(-33.80*T6**(-1.0/3.0))
    #evaluate derivative
    Texp_deriv = -2.0/3.0*T6**(-5.0/3.0)*exp(-33.80*T6**(-1.0/3.0))
    Texp_deriv = Texp_deriv+33.80/3.0*T6**(-2.0)*exp(-33.80*T6**(-1.0/3.0))
    Texp_deriv = Texp_deriv*Tsun6
    #adjust and scale
    Rd = 0.0
    Ld = 0.0
    Pd = c*rhod[2]*Texp*s
    Td = (c*rhod[3]*Texp+c*rhof*Texp_deriv)*s
    return [Rd, Ld, Pd, Td]
    
#Scaled CNO reaction
def enucCNO(P, T, M, X):
    #find rho
    rhof = rho(P, T, M, X)
    # find and interpolate mass fractions
    [fH, fHe, fOther] = X
    H = splev(M, fH)
    Other = splev(M, fOther)

    #transform T 
    T6=(T*Tsun)*1.0e-6
    #refactor and evaluate
    c = ecnoCon*H*Other
    Texp = T6**(-2.0/3.0)*exp(-152.28*T6**(-1.0/3.0))
    return c*rhof*Texp*s
    

#Scaled CNO reaction derivative  
def enucCNO_deriv(P, T, M, X):
    #find rho
    rhof = rho(P, T, M, X)
    rhod = rho_deriv(P, T, M, X)
    # find and interpolate mass fractions
    [fH, fHe, fOther] = X
    H = splev(M, fH)
    Other = splev(M, fOther)
    #transform T 
    T6=(T*Tsun)*1.0e-6
    Tsun6 = Tsun*1.0e-6
    #refactor and evaluate
    c = ecnoCon*H*Other
    Texp = T6**(-2.0/3.0)*exp(-152.28*T6**(-1.0/3.0))
    #evaluate derivative
    Texp_deriv = -2.0/3.0*T6**(-5.0/3.0)*exp(-152.28*T6**(-1.0/3.0))
    Texp_deriv = Texp_deriv+152.28/3.0*T6**(-2.0)*exp(-152.28*T6**(-1.0/3.0))
    Texp_deriv = Texp_deriv*Tsun6

    #adjust and scale
    Rd = 0.0
    Ld = 0.0
    Pd = c*rhod[2]*Texp*s
    Td = (c*rhod[3]*Texp+c*rhof*Texp_deriv)*s

    return [Rd, Ld, Pd, Td]
  

 #Scaled triple alpha reaction
def enuc3a (P, T, M, X):
     #find rho
    rhof = rho(P, T, M, X)
    # find and interpolate mass fractions
    [fH, fHe, fOther] = X
    He = splev(M, fHe)

    #transform T 
    T8=(T*Tsun)*1.0e-8
    #refactor and evaluate
    c = e3aCon*He**3
    Texp = T8**(-3.0)*exp(-44.027*T8**(-1.0))
    return c*rhof**2*Texp*s
    
#Scaled triple alpha reaction derivative  
def enuc3a_deriv(P, T, M, X):
    #find rho
    e3aCon=5.09e11
    rhof = rho(P, T, M, X)
    rhod = rho_deriv(P, T, M, X)
    # find and interpolate mass fractions
    [fH, fHe, fOther] = X
    He = splev(M, fHe)
    #transform T 

    T8=(T*Tsun)*1.0e-8
    Tsun8 = Tsun*1.0e-8

    #refactor and evaluate
    c = e3aCon*He**3
    Texp = T8**(-3.0)*exp(-44.027*T8**(-1.0))
     #evaluate derivative
    Texp_deriv = -3.0*T8**(-4.0)*exp(-44.027*T8**(-1.0))
    Texp_deriv = Texp_deriv+44.027*T8**(-5.0)*exp(-44.027*T8**(-1.0))
    Texp_deriv = Texp_deriv*Tsun8
    #adjust and scale
    Rd = 0.0
    Ld = 0.0
    Pd = c*rhod[2]*Texp*s
    Td = (c*rhod[3]*Texp+c*rhof*Texp_deriv)*s


    return [Rd, Ld, Pd, Td]
#full dL/dM (energy) equation  
def dLdM(P, T, M, X):
    #evaluate components
    enucppf = enucpp(P, T, M, X)
    enucCNOf = enucCNO(P, T, M, X)
    enuc3af = enuc3a (P, T, M, X)
    #add
    DL=enucppf+enucCNOf+enuc3af
    
    return DL

#full dL/dM (energy) derivative
def dLdM_deriv(P, T, M, X):
    #evaluate components
    enucppfd = enucpp_deriv(P, T, M, X)
    enucCNOfd = enucCNO_deriv(P, T, M, X)
    enuc3afd = enuc3a_deriv(P, T, M, X)
    #add
    Rd = enucppfd[0]+enucCNOfd[0]+enuc3afd[0]
    Ld = enucppfd[1]+enucCNOfd[1]+enuc3afd[1]
    Pd = enucppfd[2]+enucCNOfd[2]+enuc3afd[2]
    Td = enucppfd[3]+enucCNOfd[3]+enuc3afd[3]

    return [Rd, Ld, Pd, Td]    
    

