#This file calculates dT/dM according to the Scharwartzchild criterion. The related derivatives are also calculated.

#Import
from numpy import pi, minimum, maximum
from solar_kappa import opacity, opacity_deriv
from solar_P import G
from MESA import Msun, Rsun, Tsun, Lsun, Psun
    
#Constants    
k=1.3807e-16 #g.cm^2/(s^2*K)
mh=1.6733e-24 #g
c=2.99792458e10 #cm/s
a=7.566e-15 #erg.cm^-2.s^-1.K^-4=g.s^-3.K^-4 
#scaling factors 
s_conv = (Msun/Rsun**2.0)**2.0/Tsun
s2_conv = (Msun/Rsun**2.0)**2.0/Psun
s_rad = (Msun/(Rsun**2.0*Tsun**2.0))*(Lsun/(Rsun**2.0*Tsun**2.0))
    
#Calculate radiative temperature equation
def dTdM_rad(R, L, P, T, M, X):
    #find opacity
    kappaf = opacity(P, T, M, X)
    #refactor and solve
    k = 3.0/(64.0*pi**2.0*a*c)
    DTrad=-k*kappaf*L/(T**3.0*R**4.0)*s_rad
    return DTrad
#Calculate radiative temperature equation derivative
def dTraddM_deriv(R, L, P, T, M, X):
    #find opacity
    kappaf = opacity(P, T, M, X)
    kappad = opacity_deriv(P, T, M, X)
    #refactor and solve
    k = 3.0/(64.0*pi**2.0*a*c)
    #R derivative
    Rd = 4.0*k*kappaf*L/(T**3.0*R**5.0)*s_rad
    #L derivative
    Ld = -k*kappaf/(T**3.0*R**4.0)*s_rad
    #P derivative
    Pd = -k*kappad[2]*L/(T**3.0*R**4.0)*s_rad
    #T derivative
    Td = -k*kappad[3]*L/(T**3.0*R**4.0)   
    Td = Td+3.0*k*kappaf*L/(T**4.0*R**4.0)
    Td = Td*s_rad

    return [Rd, Ld, Pd, Td]
    
#Calculate convective temperature equation    
def dTdM_conv(R, P, T, M, X):


    #refactor and solve
    c = -G/(10.0*pi)
    DTconv = c*s2_conv*T*M/(P*R**4.0)
    
    return DTconv
#Calculate convective temperature equation derivative    
def dTconvdM_deriv(R, P, T, M, X): 
    #solve
    c = -G/(10.0*pi)
    #R derivative
    Rd = -4.0*(c*s2_conv)*T*M/(P*R**5)
    #L derivative
    Ld = 0.0
    #P derivative
    Pd = -1.0*(c*s2_conv)*T*M/(P**2*R**4)
    #T derivative
    Td = (c*s2_conv)*M/(P*R**4)
    
    return [Rd, Ld, Pd, Td]
#Calculate temperature equation     
def dTdM(R, L, P, T, M, X):
    #evaluate
    rad = dTdM_rad(R, L, P, T, M, X)
    conv = dTdM_conv(R, P, T, M, X)
    #Apply schwartzchild criterion 
    DT = maximum(rad, conv)

    return DT
#Calculate temperature equation  derivative
def dTdM_deriv(R, L, P, T, M, X):
    #evaluate
    conv_deriv = dTconvdM_deriv(R, P, T, M, X) #minimum(DTrad, DTconv)
    rad_deriv = dTraddM_deriv(R, L, P, T, M, X) #minimum(DTrad, DTconv)
    #evaluate
    rad = dTdM_rad(R, L, P, T, M, X)
    conv = dTdM_conv(R, P, T, M, X)
    #Apply schwartzchild criterion 
    dT_deriv = conv_deriv
    if rad > conv:
        dT_deriv = rad_deriv
                            
    return dT_deriv
