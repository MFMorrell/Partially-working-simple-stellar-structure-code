## This files extracts the values from a mesa data file and transforms them to the required form.


#Import
from pylab import *
import numpy



#Maximum Stellar values in the sun
Msun=1.989e33 #grams
Rsun=6.96e10  #cm
Lsun=3.85e33  #erg/s
Psun=2.334e17 #g/(cm s^2) (barye)
Tsun=1.573e7 #K


#Data extraction routine
def extract_MESA_values():   
    #choose file and extract
    solar=genfromtxt('solar_cal.data',skip_header=5,names=True)
     
    #Extract and scale as required (0--> 1)
    Tmesad=solar['temperature']
    Tmesa=Tmesad[::-1]
    
    Pmesad=solar['pressure']
    Pmesa=Pmesad[::-1]
    
    Mmesad=solar['mass']*Msun
    Mmesa=Mmesad[::-1]

    Rmesad=solar['radius']*Rsun
    Rmesa=Rmesad[::-1]

    Lmesad=solar['luminosity']*Lsun
    Lmesa=Lmesad[::-1]
    
    Valmesa=array([Rmesa,Lmesa,Tmesa,Pmesa])
    rhomesa = solar['logRho'][::-1]
    kappamesa=solar['log_opacity'][::-1]
    
    h = solar['h1']
    h = h[::-1]
    he =  solar['he4']
    he = he[::-1]
    gradr = solar['gradr']
    gradr = gradr[::-1]
    
    grada = solar['grada']
    grada = grada[::-1]
    
    gradT = solar['gradT']
    gradT = gradT[::-1]
    
    
    MESA_values = {}
    MESA_values['Tmesa'] = Tmesa
    MESA_values['Pmesa'] = Pmesa
    MESA_values['Mmesa'] = Mmesa
    MESA_values['Rmesa'] = Rmesa
    MESA_values['Lmesa'] = Lmesa
    MESA_values['Valmesa'] = Valmesa
    MESA_values['rhomesa'] = rhomesa
    MESA_values['H'] = h
    MESA_values['He'] = he
    MESA_values['gradr'] = gradr
    MESA_values['grada'] = grada
    MESA_values['gradT'] = gradT
    MESA_values['kappamesa'] = kappamesa
    
  
    
    
    return MESA_values

