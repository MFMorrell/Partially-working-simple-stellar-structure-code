### This files takes the equations defined in the other files and complies them into one equation. This gives us a value for the Jacobian, the equations, an initial grind, along with boundary conditions.

#import
from pylab import *
from MESA import extract_MESA_values, Lsun, Rsun, Tsun, Psun, Msun
from numpy import array
from scipy.interpolate import splrep, splev
from solar_P import dPdM, dPdM_deriv
from solar_R import dRdM, dRdM_deriv
from solar_L import dLdM, dLdM_deriv
from solar_T import dTdM, dTdM_deriv

#Define composition function
def star_composition():
    #extract values   
    H = MESA_values['H']
    He = MESA_values['He']
    Other = 1.0-H-He
    Mmesa = MESA_values['Mmesa']/Msun
    #interpolate
    fH = splrep(Mmesa, H, k = 2)
    fHe = splrep(Mmesa, He, k = 2)
    fOther = splrep(Mmesa, Other, k = 2)
   
    return [fH, fHe, fOther]

 #Define system of ODES   
def solar_function(M , x):
#readability of x
    R=x[0]    
    L=x[1]
    P=x[2]
    T=x[3]
    #get mass fractions
    X = star_composition()
    #evaluate
    Reqn = dRdM(R, P, T, M, X)
    Leqn = dLdM(P, T, M, X)
    Peqn = dPdM(R, M)
    Teqn = dTdM(R, L, P, T, M, X)
    return array([Reqn, Leqn, Peqn, Teqn])

# Define Jacobian
def solar_function_derivative(M , x):
#readability of x
    R=x[0]    
    L=x[1]
    P=x[2]
    T=x[3]
    #Mass fractions
    X = star_composition()
    #evaluate
    Rderiv = dRdM_deriv(R, P, T, M, X)
    Lderiv = dLdM_deriv(P, T, M, X)
    Pderiv = dPdM_deriv(R, M)
    Tderiv = dTdM_deriv(R, L, P, T, M, X)

    return array([Rderiv,Lderiv,Pderiv,Tderiv])
#extract values
MESA_values = extract_MESA_values()

#set range for solution in terms of indicates
def solar_set_range(start, end):
    global a
    global b
    a = start
    b = end
    
#boundary conditions are in the for array1=0 and array2=0. x0,xf are the centre and surf values.
def solar_boundary_conditions(x0,xf):
    #extract values
    Lmesa = MESA_values['Lmesa']/Lsun
    Rmesa = MESA_values['Rmesa']/Rsun
    Tmesa = MESA_values['Tmesa']/Tsun
    Mmesa = MESA_values['Mmesa']/Msun
    Pmesa = MESA_values['Pmesa']/Psun
    #set core values
    r0=Rmesa[a]
    l0=Lmesa[a] 
    #surface values
    Psurf=Pmesa[b]
    Tsurf=Tmesa[b]
    return (array([x0[0]-r0,x0[1]-l0]),  
            array([xf[2]-Psurf,xf[3]-Tsurf ]))

#define range for the mass to get initial and final values            
def solar_mass_range():
    Mmesa = MESA_values['Mmesa']/Msun
    m0=Mmesa[a]
    mtot=Mmesa[b]
    return (m0, mtot )
# generate grid for mass  
def solar_grid():
    grid = MESA_values['Mmesa']/Msun
    return grid[a:(b+1)]
#generate grid for variables    
def solar_soln_guess():
    Rmesa = MESA_values['Rmesa']/Rsun
    Lmesa = MESA_values['Lmesa']/Lsun
    Pmesa = MESA_values['Pmesa']/Psun
    Tmesa = MESA_values['Tmesa']/Tsun

    return array([Rmesa[a:(b+1)], Lmesa[a:(b+1)], Pmesa[a:(b+1)], Tmesa[a:(b+1)]])
    
