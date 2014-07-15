#This is a plotting file and contains two types of plots
#plot_solar_test4 which tests the derivatives and equations used to those of the source file
#plot_soln which compares a solution with the input data file.


from pylab import *

from MESA import extract_MESA_values
from MESA import Msun, Tsun, Rsun, Psun, Lsun
from solar_equn4 import star_composition
from solar_rho import rho
from solar_kappa import opacity
from solar_T import dTdM_rad
from solar_T import dTdM_conv
from solar_T import dTdM
from solar_P import dPdM
from solar_R import dRdM
from solar_L import dLdM



def eval_rho(MESA_values, X, rho):
    Pmesa = MESA_values['Pmesa']
    Tmesa = MESA_values['Tmesa']
    Mmesa = MESA_values['Mmesa']

    return rho(Pmesa, Tmesa, Mmesa, X)
    
def eval_kappa(MESA_values, X, kappa, rho):   
    Pmesa = MESA_values['Pmesa']
    Tmesa = MESA_values['Tmesa']
    Mmesa = MESA_values['Mmesa']
    return kappa(Pmesa, Tmesa, Mmesa, X)
    




def plot_rho(i, MESA_values, X, rho):
    plt.figure(i)
    Mmesa = MESA_values['Mmesa']
    rhomesa = MESA_values['rhomesa']
    an_rho = eval_rho(MESA_values, X, rho)
    plot1=plt.plot(log10(Mmesa),log10(an_rho))
    plot2=plt.plot(log10(Mmesa),rhomesa,'+')
    plt.title("Density vs Mass")
    plt.ylabel("Log10(Density -- rho )")
    plt.xlabel("Log10(Mass)")
    plt.legend((plot1[0],plot2[0]),('Analytical', 'MESA'))
    plt.figure(i+1)
    delrho=abs(10**rhomesa-an_rho)/abs(10**rhomesa)
    plt.plot(log10(Mmesa),delrho)
    plt.title("Relative Error in Density")
    plt.ylabel("Relative Error in Density -- rho")
    plt.xlabel("Log10(Mass)")
    
def plot_kappa(i, MESA_values, X, kap, rho):
    plt.figure(i)
    Mmesa = MESA_values['Mmesa']
    Tmesa = MESA_values['Tmesa']
    kappamesa= MESA_values['kappamesa']
    an_kappa = eval_kappa(MESA_values, X, kap, rho)
    plot1=plt.plot(log10(Tmesa),log10(an_kappa))
    plot2=plt.plot(log10(Tmesa), kappamesa,'+')
    plt.title("Opacity vs Tempr")
    plt.legend((plot1[0],plot2[0]),('Analytical', 'MESA'))
    plt.ylabel("Log10(Opacity -- kappa )") #cm^2/g
    plt.xlabel("Log10(Temperature)")
    plt.figure(i+1)
    delkappa=abs(10**kappamesa-an_kappa)/abs(10**kappamesa)
    plt.plot(log10(Tmesa),delkappa)
    plt.title("Relative Error in Opacity")
    plt.ylabel("Error in Opacity -- kappa")
    plt.xlabel("Log10(Temperature)")
    


def plot_dRdM(i, MESA_values, X, DRDM):
    plt.figure(i)
    Mmesa = MESA_values['Mmesa']
    Rmesa = MESA_values['Rmesa']
    Pmesa = MESA_values['Pmesa']
    Tmesa = MESA_values['Tmesa']
    
    
    NumdR=np.diff(Rmesa)/np.diff(Mmesa)
    an_DR = DRDM(Rmesa, Pmesa, Tmesa, Mmesa, X)
    plot1=plt.plot(log10(Mmesa),log10(an_DR))
    plot2=plt.plot(log10(Mmesa[0:-1]),log10(NumdR),'+')
    plt.title("Radius Gradient vs Mass")
    plt.ylabel("Log10(Radius Gradient -- DR/DM )")#cm/g
    plt.xlabel("Log10(Mass)")
    plt.legend((plot1[0],plot2[0]),('Analytical', 'MESA'))
    plt.figure(i+1)
    delDR=abs(NumdR-an_DR[0:-1])/abs(NumdR)
    plt.plot(log10(Mmesa)[0:-1],delDR)
    plt.title("Relative Error in Radius Gradient")
    plt.ylabel("Relative Error in Radius Gradient -- DR/DM")
    plt.xlabel("Log10(Mass)")
    
    
    
    
def plot_dLdM(i, MESA_values, X, DLDM):
    plt.figure(i)
    Pmesa = MESA_values['Pmesa']
    Lmesa = MESA_values['Lmesa']
    Tmesa = MESA_values['Tmesa']
    Mmesa = MESA_values['Mmesa']
    NumdL=np.diff(Lmesa)/np.diff(Mmesa)
    an_DL = DLDM(Pmesa, Tmesa, Mmesa, X)
    plot1=plt.plot(log10(Mmesa),log10(an_DL))
    plot2=plt.plot(log10(Mmesa[0:-1]),log10(NumdL),'+')
    ylim([-10,10]) 
    plt.title("Energy vs Mass")
    plt.ylabel("Log10(Energy -- DL/DM )")
    plt.xlabel("Log10(Mass)")
    plt.legend((plot1[0],plot2[0]),('Analytical', 'MESA'))
    plt.figure(i+1)
    
    delDL=abs(NumdL-an_DL[0:-1])/abs(NumdL)
    plt.plot(log10(Mmesa)[0:-1],delDL)
    ylim([0,1.1]) 
    plt.title("Relative Error in Energy")
    plt.ylabel("Relative Error in Energy -- dL/dM")
    plt.xlabel("Log10(Mass)")
    


def plot_dPdM(i, MESA_values, DPDM):
    plt.figure(i)
    Pmesa = MESA_values['Pmesa']
    Mmesa = MESA_values['Mmesa']
    Rmesa = MESA_values['Rmesa']
    

    NumdP=np.diff(Pmesa)/np.diff(Mmesa)
    an_DP = DPDM(Rmesa, Mmesa)
    plot1=plt.plot(log10(Mmesa),log10(-an_DP))
    plot2=plt.plot(log10(Mmesa[0:-1]),log10(-NumdP),'+')
    plt.title("Pressure Gradient vs Mass")
    plt.ylabel("Log10(Pressure Gradient -- dP/dM )") #Ba/g
    plt.xlabel("Log10(Mass)")
    plt.legend((plot1[0],plot2[0]),('Analytical', 'MESA'))
    plt.figure(i+1)
    delDP=abs(NumdP-an_DP[0:-1])/abs(NumdP)
    plt.plot(log10(Mmesa)[0:-1],delDP)
    plt.title("Relative Error in Pressure Gradient")
    plt.ylabel("Relative Error in Pressure Gradient --  dP/dM ")
    plt.xlabel("Log10(Mass)")
    


    
def plot_dTdM(i, MESA_values, X, DTDM):
    Lmesa = MESA_values['Lmesa']
    Rmesa = MESA_values['Rmesa']
    Tmesa = MESA_values['Tmesa']
    Pmesa = MESA_values['Pmesa']
    Mmesa = MESA_values['Mmesa']

    plt.figure(i)
    NumdT=np.diff(Tmesa)/np.diff(Mmesa)
    an_DT = DTDM(Rmesa, Lmesa, Pmesa, Tmesa, Mmesa, X)
    plot1=plt.plot(log10(Mmesa),log10(-an_DT))
    plot2=plt.plot(log10(Mmesa[0:-1]),log10(-NumdT),'+')
    plt.title("Temperature  Gradient vs Mass")
    plt.legend((plot1[0],plot2[0]),('Analytical', 'MESA'))
    plt.ylabel("Log10(Temperature Gradient -- dT/dM)") #K/g
    plt.xlabel("Log10(Mass)")
    plt.figure(i+1)
    delDT=abs(NumdT-an_DT[0:-1])/abs(NumdT)
    plt.plot(log10(Mmesa)[0:-1],delDT)
    plt.title("Relative Error in Temperature Gradient")
    plt.ylabel("Relative Error in Temperature Gradient -- dT/dM ")
    plt.xlabel("Log10(Mass)")
 

def plot_dTdMcomp(i, MESA_values,X, DTradDM,DTconvDM):
    plt.figure(i)
    Lmesa = MESA_values['Lmesa']
    Pmesa = MESA_values['Pmesa']
    Rmesa = MESA_values['Rmesa']
    Tmesa = MESA_values['Tmesa']
    Mmesa = MESA_values['Mmesa']
    NumdT=np.diff(Tmesa)/np.diff(Mmesa)
    an_DTrad = DTradDM(Rmesa, Lmesa, Pmesa, Tmesa, Mmesa, X)
    an_DTconv = DTconvDM(Rmesa,Pmesa,Tmesa, Mmesa, X)
    plot1=plt.plot(log10(Mmesa),log10(-an_DTconv))
    plot2=plt.plot(log10(Mmesa),log10(-an_DTrad),'*')
    plot3=plt.plot(log10(Mmesa[0:-1]),log10(-NumdT),'+')
    plt.legend((plot1[0],plot2[0],plot3[0]),('Analytical Convective Gradient','Analytical Radiative Gradient', 'MESA Full Gradient'))
    plt.title("Temperature Gradient Type Comparison")
    plt.ylabel("Log10(Temperature Gradient -- dT/dM)")
    plt.xlabel("Log10(Mass)")


    
#Test the derivatives and equations used to those of the source file
    
def plot_solar_test4():
    
    MESA_values= extract_MESA_values()
    MESA_values['Mmesa']  =  MESA_values['Mmesa'] /Msun
    MESA_values['Tmesa']  =  MESA_values['Tmesa'] /Tsun
    MESA_values['Rmesa']  =  MESA_values['Rmesa'] /Rsun
    MESA_values['Pmesa']  =  MESA_values['Pmesa'] /Psun
    MESA_values['Lmesa']  =  MESA_values['Lmesa'] /Lsun
    X = star_composition()
   
    plot_dRdM(1, MESA_values, X, dRdM)
    plot_dLdM(3, MESA_values, X, dLdM)
    plot_dPdM(5, MESA_values, dPdM)    
    plot_dTdM(7, MESA_values, X, dTdM)
    plot_dTdMcomp(9, MESA_values, X, dTdM_rad, dTdM_conv)# 1 plot
    plot_rho(10, MESA_values, X, rho)
    plot_kappa(12, MESA_values, X, opacity, rho)


    show()

#Compare solution with the input data file.
def plot_soln(M, R, L, P, T, a, b):
    MESA_values= extract_MESA_values()
    Mmesa = MESA_values['Mmesa'][a:b]/Msun 
    plt.figure(1)     
    Rmesa = MESA_values['Rmesa'][a:b]/Rsun
    plot1=plt.plot(Mmesa, Rmesa, '+')
    plot2=plt.plot(M, R)
    plt.legend((plot1[0],plot2[0]),('MESA','Result' ))
    plt.title("Radius Comparison")
    plt.ylabel("Radius")
    plt.xlabel("Mass")
    
    plt.figure(2)
    Lmesa = MESA_values['Lmesa'][a:b]/Lsun
    plt.plot(Mmesa, Lmesa,'+')
    plt.plot(M, L)
    plt.legend((plot1[0],plot2[0]),('MESA','Result' ))
    plt.title("Luminosity Comparison")
    plt.ylabel("Luminosity")
    plt.xlabel("Mass")


    plt.figure(3)
    Pmesa = MESA_values['Pmesa'][a:b]/Psun
    plt.plot(Mmesa, Pmesa,'+')
    plt.figure(3)
    plt.plot(M, P)
    plt.legend((plot1[0],plot2[0]),('MESA','Result' ))
    plt.title("Pressure Comparison")
    plt.ylabel("Pressure")
    plt.xlabel("Mass")

    plt.figure(4)
    Tmesa = MESA_values['Tmesa'][a:b]/Tsun
    plt.plot(Mmesa, Tmesa,'+')
    plt.plot(M, T)
    plt.legend((plot1[0],plot2[0]),('MESA','Result' ))
    plt.title("Temperature Comparison")
    plt.ylabel("Temperature")
    plt.xlabel("Mass")
    show()
    

