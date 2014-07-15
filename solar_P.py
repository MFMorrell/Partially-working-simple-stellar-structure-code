#This file calculates dP/dM and its associated derivative

#Import
from numpy import pi, exp
from MESA import Psun, Msun, Rsun
#Constant
G=6.674e-8 # units of cm^3/(g.s^2)
#scaling factor
s = (Msun/Rsun**2.0)**2.0/Psun

#Pressure
#Unscaled version
def dPdM_orig(R, M):
    DP=-G*M/(4.0*pi*(R**4))
    
    return DP
#Unscaled version    
def dPdM_deriv_orig(R, M):
 
   
    Rd = 4.0*G*M/(4.0*pi*(R**5))
    Ld = 0.0
    Pd = 0.0
    Td = 0.0
    
    return [Rd, Ld, Pd, Td]
 #dP/dM equation  
def dPdM(R, M):
    
    DP=(-G*s)*M/(4.0*pi*(R**4.0))
    
    return DP
 #derivatives of dP/dM for the Jacobian    
def dPdM_deriv(R, M):
 
    #R derivative   
    Rd = 4.0*G*M/(4.0*pi*(R**5.0))*s
    #L derivative
    Ld = 0.0
    #P derivative
    Pd = 0.0
    #T derivative
    Td = 0.0
        
    
    return [Rd, Ld, Pd, Td]
