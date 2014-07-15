# This file contains the calculations for opacity and it derivative. It also contains the code to extract from opal. Please note that opal MUST be turned into a python module by running
#>f2py -c -m opal xztrin21.f
# before running this file.

#Required imports
from solar_rho import rho, rho_deriv
from MESA import Tsun
from scipy.interpolate import splrep, splev
from numpy import array
#Old constants
kffconst=3.68e22 #cm^5*k^3.5/g^2
kbfconst=4.34e25# cm^5*k^3.5/g^2      "
kesconst=0.2 #cm^2/g

s = Tsun**-3.5 

def opacity_table(T, X, rhof):
    #import opal calculations
    from opal import  opacgn93, e
    
    [H, He, Other] = X

#    print opacgn93.__doc__
#    print e.__doc__
#Scale
    T6 = T*(Tsun*1.0E-6)
    S = rhof/T6**3
    
 #   if (R < 2.4E-9):
 #       kap = kesconst*(1.0+H)
 #   else:
 #Run Opal
    opacgn93(Other, H, T6, S)
#Convert output    
    kap = 10.0**e.opact

    return kap
    
    
#calculate opacity    
def opacity(P, T, M, X):   
    #import mass fractions and set to mass grid  
    [fH, fHe, fOther] = X
    H = splev(M, fH)
    He = splev(M, fHe)
    Other = splev(M, fOther)
    #find rho
    rhof = rho(P, T, M, X)
    # error handling
    if isinstance(T, float):
        kap = opacity_table(T, [H, He, Other], rhof)
    else:
        n = len(T)
        kap = [0]*n
        for i in range(0, n):
            kap[i] = opacity_table(T[i], [H[i], He[i], Other[i]], rhof[i])
        kap = array(kap)

 #   kapES=kesconst*(1.0+H)
    
 #   kapFF=kffconst*(He+H)*(1.0+H)*rhof*T**-3.5*s 
 #   kapBF=kbfconst*(Other)*(1.0+H)*rhof*T**-3.5*s
 #   kap=kapES+kapBF+kapFF 
    
    return kap
#approximate derivate. Possibly slightly incorrect.
def opacity_deriv(P, T, M, X):

    h = 1.0E-8
    
    k = opacity(P, T, M, X)
    kp = opacity(P+h, T, M, X)
    kT = opacity(P, T+h, M, X)
    Rd = 0.0 
    Ld = 0.0
    Pd = (kp-k)/h
    Td = (kT-k)/h
        
    return [Rd, Ld, Pd, Td]
 

