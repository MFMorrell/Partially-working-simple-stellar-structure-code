#This file calculates dR/dm and its derivatives.

#Import 
from numpy import pi, exp
from solar_rho import rho, rho_deriv
from MESA import Msun, Rsun
#Scaling factor
s = Msun/Rsun**3.0

#dR/dM evaluation    
def dRdM(R, P, T, M, X):
 #calculate rho
    rhof = rho(P, T, M, X)
#factor equation    
    c = 4.0*pi
#evaluate. Note s is a scaling factor 
    dR = 1.0/(c*R**2.0*rhof)*s
    
    return dR
#dR/dM row of the Jacobian    
def dRdM_deriv(R, P, T, M, X):
    #calculate rho and its derivative
    rhof = rho(P, T, M, X)
    rho_derivs = rho_deriv(P, T, M, X)
    #factor equation 
    c = 4.0*pi
    #R derivative
    Rd = -2.0/(c*R**3.0*rhof)*s
    #L derivative
    Ld = 0.0
    #P derivative
    Pd = -1.0/(c*R**2.0*rhof**2.0)*rho_derivs[2]*s
    #T derivative
    Td = -1.0/(c*R**2.0*rhof**2.0)*rho_derivs[3]*s
    
    
    return [Rd, Ld, Pd, Td]
