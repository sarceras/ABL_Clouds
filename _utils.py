import numpy as np
from _constants import *

# Functions for Carbon radiative forcing.
def R_c(DeltaC): # input in tonneCarbon tC, come Kirschnbaum
    """Radiative forcing for an DeltaC"""
    return 5.35*(np.log(1 + DeltaC/(C0*2.123*10**9)))*510*10**(12)#*86400*365 # [W/m^2/year]  
#per year/tC whole earth, instantaneous value on earth

def gaCalc(h):
    """Calculate atmospheric conductance (m/s) from plant height (m)"""
    return kvc**2.*u_wind/np.log(((h+1.) - 0.64*h)/(0.13*h))**2. # change z_wind with h+1

def RF_to_Ceq(RF): # input in W/m^2
    """Translate radiative forcing into carbon stock equivalent"""
    return 2.123*10**9*C0*(np.exp(RF/(5.35*512*10**12)) - 1) # output in tC/m^2

def fDelta(T):
    Delta=3.31863*10**6*np.exp(5425.16*(0.00366099 - 1./T) )/T**2
    return Delta

def fes(T):
    es = 611.71*np.exp(2.501/0.461*10**3* (1./273.15-1./T) )
    return es

def JarvisfPhi(Q):
    fPhi = 1. - np.exp(-0.005*Q)
    return fPhi

def JarvisfPsi_l(Psi_l):
    if Psi_l < -4.5:
        fPsi_l = 0.
    elif Psi_l < -0.05:
        fPsi_l = (Psi_l + 4.5)/(-0.05 + 4.5)
    else:
        fPsi_l = 1.
    return fPsi_l

def JarvisfD(VPD):
    fD = 1./(1 + VPD/1250)
    return fD

def JarvisfTa(Ta):
    fTa = 1.-0.0016*(Ta - 298)**2
    return fTa

def gaCalc(h):
    """calculate atmospheric conductance (m/s) from plant height (m)"""
    return kvc**2.*u_wind/np.log(((h+1.)-0.64*h)/(0.13*h))**2.     #change z_wind with h+1

def ArdenBuck(T):
    es= 6.1121 * np.exp( (18.678-T/234.5) * (T/(257.14+T))) #hPa (saturation vapor pressure)
    es = es*100 #Pa
    return es   
























