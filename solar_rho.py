#The calculates the density of the star and its derivative. We note that this assume the area is fully ionized

#import
from MESA import Psun, Tsun
from scipy.interpolate import splrep, splev

#weights
mmH=1.0 #g
mmHe=4.0 #g
mmOth=12.0 #g
   
#Constants 
Rconst=8.3144621e7#cm^2.g/(K.s^2.mol)
#scaling factor
s = Psun/Tsun
    
#calculate density    
def mu(M, X):
    #get values
    [fH, fHe, fOther] = X
    #fix to mass coordinate
    H = splev(M, fH)
    He = splev(M, fHe)
    Other = splev(M, fOther)
    #Calculate
    mu = 4.0/(2.0+6.0*H+He)
    return mu
 #calculate rho      
def rho(P, T, M, X):
    #factorized to avoid round-off
    c = mu(M, X)/Rconst
    #s is scaling factor
    rhof = c*P/T*s
    
    return rhof   
 
#derivative of rho for the analytical Jacobian 
def rho_deriv(P, T, M, X):
    
    c = mu(M, X)/Rconst
    
    Rd = 0.0
    Ld = 0.0
    Pd = c/T*s
    Td = -1.0*c*P/T**2*s
        
    return [Rd, Ld, Pd, Td]
