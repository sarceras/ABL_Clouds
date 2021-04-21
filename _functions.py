import numpy as np
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy import integrate
import pandas as pd
from _constants import *
from _utils import *

def Evap(theta, q, Q, s, specie, soil):
   
    if Q <= 20:     
        Q = 20.

    lambda_ = 2.45*10**6
    
    try:  
        Psi_l = fsolve(f_solvePsil, -1, args = (theta, q, Q, s, specie, soil)) 
    except:
        raise ValueError("Error.")   

    Eva, _, _ = f_solvePsilEvap(Psi_l, theta, q, Q, s, specie, soil)
    E = Eva *rhow # kg/m^2/s
    
    H = Q - lambda_*E # Now they are both in J/m^2/s or W/m^2

    if H < 10:    
        H = 10.;
        E = (Q-H)/lambda_
    return E    

def f_solvePsil(x, theta, q, Q, s, specie, soil):
    # Define 
    d=2.
    c=2.
    LAI_=LAI[specie]
    RAIstar=RAIW[specie] #max per well watered
    Zr_=Zr[specie]  #m
    gpmax_=gpmax[specie] #um/s
    #h=height[specie]
    #gba=gaCalc(h)  #m/s   
    gba=ga[specie]
    gsmax=0.025   #m/s
    a=8.
    b_=b[soil]
    mPsi_s=mpsi_s[soil]
    Ks_=Ks[soil]/24./3600./100.  #to have it in m/s
    #print(Ks_,gba, LAI_, Zr_,  )
    gamma=cp*P0/0.622/lambda_
    VPD=fes(theta)-q*P0/0.622;
    gs=gsmax*JarvisfPhi(Q)*JarvisfPsi_l(x)*JarvisfD(VPD)*JarvisfTa(theta)
    Kroot=Ks_*s**(2*b_+3)
    RAI=RAIstar*s**(-a)   
    gsr=Kroot*RAI**0.5/g/np.pi/rhow/Zr_*10**6
    gp=gpmax_*( 1+np.exp(-(-x/d)**c) ) * 10**(-6)
    gsrp=LAI_*gsr*gp/(gsr+LAI_*gp)
    epsilon_r=fDelta(theta)
    
    Etop=(epsilon_r*Q+rho*cp*gba*VPD)*gs*LAI_;
    Ebottom=rhow*lambda_*( gamma*(gba+gs*LAI_)+gs*LAI_*epsilon_r );
    E=Etop/Ebottom
    
    Psi_s=mPsi_s*s**(-b_)
    E2=gsrp*(Psi_s-x)
    y=E-E2
    #print('y=',y)
    return y

def f_solvePsilEvap(Psi_l, theta, q, Q, s, specie, soil):
    d=2.
    c=2.
    LAI_=LAI[specie]
    RAIstar=RAIW[specie]
    Zr_=Zr[specie]  #m
    gpmax_=gpmax[specie] #um/s
    #h=height[specie]
    #gba=gaCalc(h)  #m/s   #m/s  
    gba=ga[specie]
    gsmax=0.025   #m/s
    a=8.
    b_=b[soil]
    mPsi_s=mpsi_s[soil]
    Ks_=Ks[soil]/24./3600./100.  #to have it in m/s
    
    VPD=fes(theta)-q*P0/0.622;
    gs=gsmax*JarvisfPhi(Q)*JarvisfPsi_l(Psi_l)*JarvisfD(VPD)*JarvisfTa(theta)
    #print('gs=', gs)
    Kroot=Ks_*s**(2*b_+3.)
    RAI=RAIstar*s**(-a)   
    gsr=Kroot*RAI**0.5/g/np.pi/rhow/Zr_*10**6
    #print('gsr=', gsr)
    gp=gpmax_*(1+np.exp(-(-Psi_l/d)**c) ) * 10**(-6)  #m/MPa/s
    gsrp=LAI_*gsr*gp/(gsr+LAI_*gp)   #m/s
    Psi_s=mPsi_s*s**(-b_)
    Eva=gsrp*(Psi_s-Psi_l) #m/s
    return Eva, gs, VPD

def simpleLCL( T0,q0 ):
    Taud=-g/cp
    epsilon=Rd/Rv
    
    TT=np.arange(230,351,1)
    es=ArdenBuck(TT-273.15)
    qs=epsilon*es/P0
    Tdew=interp1d(qs,TT)(q0)
    zLCL=125*(T0-Tdew)
    TLCL=T0+zLCL*Taud
    
    return zLCL, TLCL

def thetav2TP(FA_params):
    #transfer theta profile to T and P profile

    z = np.arange(0, 18010, 10) # for calculating CAPE

    q = FA_params["yq"]*z + FA_params["q0"]
    q[q<=0] = 0

    P = np.empty(z.shape)
    T = np.empty(z.shape)
    Thetav = FA_params["ythv"]*z + FA_params["thv0"]
    Theta = Thetav/(1 + 0.61*q)

    P[0] = P0
    T[0] = Theta[0]
    dz = z[1] - z[0]
    
    for i in range(len(P)-1):
        T[i] = Theta[i]*(P[i]/P0)**(Rd/cp)
        dPdz = -P[i]*g/Rd/T[i]
        P[i+1] = P[i] + dPdz*dz
    
    T[-1] = Theta[-1]*(P[-1]/P0)**(Rd/cp)
    
    return z, P, T, q

def Adiabat( TLCL,PLCL ):

    epsilon=Rd/Rv

    #moist adiabatic
    Pm=np.linspace(PLCL,20000,500)
    Tm=np.empty(Pm.shape)
    Tm[0]=TLCL
    for i in range(len(Pm)-1):
        es=ArdenBuck(Tm[i]-273.15)
        qs=epsilon*es/Pm[i]
        dTdP1=Rd*Tm[i]+lambda_*qs;
        dTdP2=(  cp+lambda_**2*qs/Rv/Tm[i]/Tm[i]  ) * Pm[i]
        dTdP=dTdP1/dTdP2
        Tm[i+1]=Tm[i]+dTdP* (Pm[i+1]-Pm[i])

    #dry adiabatic
    Pd=np.linspace(P0,PLCL-1,100)
    T0=TLCL*(P0/PLCL)**(Rd/cp)
    Td=T0*(Pd/P0)**(Rd/cp)
    
    #combined
    PAD=np.concatenate((Pd,Pm))
    TAD=np.concatenate((Td,Tm))
    
    return TAD, PAD

def getAllV(Tsrd,Psrd,TAD,PAD,qsrd,q,PLCL):
    Rd = 287.
    #synchronize the surrouding index
    Tsrd=interp1d(Psrd,Tsrd)(PAD)
    qsrd=interp1d(Psrd,qsrd)(PAD)
    Psrd=PAD

    # virtual potential temperature
    index1=(PAD>PLCL) # below LCL
    index2=(PAD<=PLCL)  # above LCL
    es=ArdenBuck(TAD[index2]-273.15)
    qs=0.622*es/PAD[index2]
    TADV=np.zeros(index2.shape)
    TADV[index2]=TAD[index2]*(1+0.61*qs)
    TADV[index1]=TAD[index1]*(1+0.61*q)
    
    TsrdV=Tsrd*(1+0.61*qsrd)

    # find climate index
    A=TADV-TsrdV

    index=np.logical_and(A>0, PAD<PLCL)
    if sum(index)>2:
        X=A[index]
        PP=PAD[index]
        CAPE=-np.trapz(X*Rd,np.log(PP))
        LFC=PP[0]
        LNB=PP[-1]

        index=(PAD>LFC)
        X=A[index]
        PP=PAD[index]
        if len(PP)<=1:
            CIN=0
        else:
            CIN=np.trapz(X*Rd,np.log(PP));
        
    else:
        
        CAPE=np.nan
        LFC=np.nan
        LNB=np.nan
        CIN=np.nan
    
    return CAPE,CIN, LFC,LNB    

def odefunV(t, y, FA_params, QQ, tt, s, specie, soil):
    # Define matrix for solutions of ode.
    dy = np.zeros((len(y)))
    
    h = y[0]
    thetav = y[1]
    q = y[2]
    
    theta = thetav/(1. + 0.61*q)
    
    Q = interp1d(tt, QQ, fill_value = "extrapolate")(t) # Extrapolate takes value even out the time range.
    
    if Q <= 20:     
        Q = 20.
  
    try:  
        Psi_l = fsolve(f_solvePsil, x0 = -1, args = (theta, q, Q, s, specie, soil)) 
    except:
        raise ValueError("Error.")   

    # Solve energy balance.
    Eva, _, _ = f_solvePsilEvap(Psi_l, theta, q, Q, s, specie, soil)
    
    E = Eva*rhow # kg/m^2/s
    H = Q-lambda_*E # W/m^2
    le = lambda_*E
    
    if H < 10:    
        H = 10.;
        E = (Q-H)/lambda_
        le = lambda_*E

    Hv = H + 0.61*theta*cp*E
    
    dy[0] = (1 + 2*beta)*Hv/rho/cp/FA_params["ythv"]/h
    dy[1] = (Hv + rho*cp*(FA_params["thv0"] + FA_params["ythv"]*h - thetav)*dy[0]) / (rho*cp*h)
    dy[2] = (E + rho*(FA_params["q0"] + FA_params["yq"]*h - q)*dy[0]) / (rho*h)

    return dy

def AP(FA_params, Q_max, s, specie, soil, alfa, tt):
    # Define list of times from sunrise to sunset.
    tt = tt*3600 # s
    # Define sunrise.
    t0 = tt[0] # s
    # Define parabola radiation.
    QQ = Q_max*(1 - ((tt - t0)/t0 - 1)**2)

    QQ = QQ*(1.-alfa)

    # Define initial conditions.
    # Height ABL.
    h0 = 50. # m
    # FA params.
    thetav0 = FA_params["ythv"]*(h0/2) + FA_params["thv0"]
    q0 = FA_params["yq"]*(h0/2) + FA_params["q0"]
    
    y0 = [h0, thetav0, q0]
    
    # Set for integration.
    tspan = (tt[0], tt[-1])

    sol = integrate.solve_ivp(lambda t, y: odefunV(t, y, FA_params, QQ, tt, s, specie, soil), 
                              t_span = tspan, y0 = y0, t_eval = tt)

    h = sol.y[0]
    thetav = sol.y[1]
    q = sol.y[2]
    T=sol.t 

    theta=thetav/(1.+0.61*q) 

    CAPE=np.empty(theta.shape)
    CIN=np.empty(theta.shape)
    LFC=np.empty(theta.shape)
    LNB=np.empty(theta.shape)
    zLCL=np.empty(theta.shape)
    TLCL=np.empty(theta.shape)
    
    #define lists to store flux values for diurnal evolution
    LE=np.empty(theta.shape)
    H=np.empty(theta.shape)
    
    for i in range(len(theta)):
        lE=float(lambda_*Evap(theta[i],q[i],QQ[i],s, specie, soil))
        LE[i]=lE
        H[i]=QQ[i]-lE

        _, Psrd, Tsrd,qsrd = thetav2TP(FA_params) #Psrd independent of y_q
   
        try:
            zLCL, TLCL = simpleLCL( theta,q ) 
            PLCL=P0*(TLCL/theta)**(cp/Rd)
            Ph=P0*((theta-g/cp*h)/theta)**(cp/Rd)
            # adiabatic lifting
            TAD, PAD = Adiabat( TLCL[i],PLCL[i] )

            Tsrd=interp1d(Psrd,Tsrd)(PAD)  
            qsrd=interp1d(Psrd,qsrd)(PAD)
            Psrd=PAD
            # synchronize the surrouding below ABL
            index=Psrd>Ph[i]
            Tsrd[index]=TAD[index]
            qsrd[index]=q[i]

            CAPE[i],CIN[i], LFC[i],LNB[i]  = getAllV(Tsrd,Psrd,TAD,PAD,qsrd,q[i],PLCL[i])
        except:
            pass

    # transfer pressure altitude
    try:
        zsrd, Psrd, _,_ = thetav2TP( gamma_thetav, thetavf0, gamma_q,qf0 )
        LFC=interp1d(Psrd,zsrd)(LFC)
        LNB=interp1d(Psrd,zsrd)(LNB)
    except:
        pass

    return T, h, zLCL, CAPE, CIN, LFC, LNB, theta, q, LE,H, QQ, # alfa

def daily_sim(FA_params, Q_max, s, cloud_albedo, specie, soil, alfa, tt):    
    # Compute AP.
    T, h, zLCL, CAPE, CIN, LFC, LNB, theta, q, LE, H, QQ  = AP(FA_params, Q_max, s, specie, soil, alfa, tt)

    interp0 = interp1d(T, QQ, fill_value = "extrapolate")
    T0, T1 = fsolve(interp0, [20000, 70000])
    t = np.linspace(T0, T1) 

    area_Q = np.trapz(interp0(t), t)
    mean_area_Q = area_Q/(24.*3600.)
    
    interp1 = interp1d(T, h)
    interp2 = interp1d(T, zLCL)
    
    def difference(x):
        return np.abs(interp1(x) - interp2(x))

    try:
        t_at_crossing = float(fsolve(difference, 30000))
    except:   
        t_at_crossing = np.nan

    if np.isnan(t_at_crossing): 
        y_at_crossing = np.nan 
        Cape_cross = np.nan
        area_Q_cloud = area_Q 
    else:
        y_at_crossing = interp1(t_at_crossing)
        Cape_cross = interp1d(T, CAPE)(t_at_crossing)
        x1 = np.linspace(T0, t_at_crossing)  
        x2 = np.linspace(t_at_crossing, T1)
        area_Q_cloud = np.trapz(interp0(x1), x1) + np.trapz(interp1d(T, (1.-cloud_albedo)*QQ, 
                                                                     fill_value='extrapolate')(x2), x2) 

    mean_area_Q_cloud = area_Q_cloud/(24.*3600.)
    area_LE = np.trapz(interp1d(T, LE)(t), t)
    mean_area_LE = area_LE/(24.*3600.)

    #just for daily graphs with clouds
    QQ_cloud = []
    for i in range(len(T)):
        if np.isnan(t_at_crossing):
            old_Q = QQ[i]
            new_Q = old_Q
        elif T[i] <= t_at_crossing:
            old_Q = QQ[i]
            new_Q = old_Q
        else:
            old_Q = QQ[i]
            new_Q = old_Q*(1.-cloud_albedo)
        QQ_cloud.append(new_Q)  

    QQ_cloud = np.array(QQ_cloud)  

    return pd.Series([T, QQ, QQ_cloud, t_at_crossing, Cape_cross, mean_area_Q, mean_area_Q_cloud, mean_area_LE], 
                     index = ["T", "Q", "Q_cloud", "Crossing_time", "CAPE_at_crossing", "mean_daily_Q", "mean_daily_Q_cloud", 
                              "mean_daily_LE"], name = "Daily_param") 