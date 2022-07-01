# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 11:47:24 2022

@author: hidde
"""



import numpy as np
import matplotlib.pyplot as plt

#starting april 25th
deathsperweek=[2,0,9,3,14,17,43,112,168,267,470,725,1089,1843,2010,2817,3880,4237,6102,6988,6544,7165,5533,4929,4327,2665,1421,1031,1414,1050,652,333,210,243,281]
deathsperweek=[0]*40+deathsperweek
totaldeaths=np.cumsum(deathsperweek)

#end december 19th

#plt.plot(deathsperweek)
tlst=(np.linspace(0,74,75)*7)

population=460000
h=0.001
T=2
N=int(T/h)

def loop(rho):
    Sh=np.zeros(N)
    Ih=np.zeros(N)
    Rh=np.zeros(N)
    Dh=np.zeros(N)
    
    F=np.zeros(N)
    Nf=np.zeros(N)
    
    Sr=np.zeros(N)
    Ir=np.zeros(N)
    Rr=np.zeros(N)
    Dr=np.zeros(N)
    Vr=np.zeros(N)
    
    Tr=np.zeros(N)
    
    s=np.zeros(N)
    t=np.zeros(N)
    
    Sh[0]=population
    Sr[0]=55000
    Rr[0]=30000
    Ir[0]=500
    Nf[0]=5
    
    
    r_r=5
    p=0.975
    K_r=0.2*population
    #0.5<aK_r<20
    a=10/K_r
    d_r=0.2
    beta_r=4.7
    gamma_r=20
    g_r=0.02
    r_f=20
    d_f=10
    K_f=6.57
    mu_f=10
    beta_h=0.01
    gamma_h=26
    mu_h=0.5
    
    for i in range (0,round(7*40*1000/365),1):
        Tr[i]=Sr[i]+Ir[i]+Rr[i]+Vr[i]
        s[i]=np.exp(-a*Tr[i])
        t[i+1]=t[i]+h
    
        Sr[i+1]=Sr[i]+h*(r_r*Sr[i]*(1-Tr[i]/K_r)+r_r*Rr[i]*(1-p)-beta_r*Sr[i]/(Tr[i]) *F[i]*(1-s[i])-d_r*Sr[i])
        Ir[i+1]=Ir[i]+h*(beta_r*Sr[i]/Tr[i]*F[i]*(1-s[i])-(gamma_r+d_r)*Ir[i])
        Rr[i+1]=Rr[i]+h*(r_r*Rr[i]*(p-Tr[i]/K_r)+gamma_r*g_r*Ir[i]-d_r*Rr[i])
        Dr[i+1]=Dr[i]+h*(+gamma_r*(1-g_r)*Ir[i]+d_r*Tr[i])
    
        
        Nf[i+1]=Nf[i]+h*(r_f*Nf[i]*(1-Nf[i]/K_f)+mu_f/Tr[i]*F[i]*(1-s[i]))
        F[i+1]=F[i]+h*((d_r+gamma_r*(1-g_r))*Ir[i]*Nf[i]-d_f*F[i])
        
        Sh[i+1]=Sh[i]+h*(-beta_h*Sh[i]*F[i]*s[i])
        Ih[i+1]=Ih[i]+h*(beta_h*Sh[i]*F[i]*s[i]-gamma_h*Ih[i])
        Rh[i+1]=Rh[i]+h*(gamma_h*Ih[i]*(1-mu_h))
        Dh[i+1]=Dh[i]+h*(gamma_h*Ih[i]*mu_h)
    
    for i in range (round(7*40*1000/365),N-1,1):
        Tr[i]=Sr[i]+Ir[i]+Rr[i]+Vr[i]
        s[i]=np.exp(-a*Tr[i])
        t[i+1]=t[i]+h
        Sr[i+1]=Sr[i]+h*(r_r*(Vr[i]+Sr[i])*(1-Tr[i]/K_r)+r_r*Rr[i]*(1-p)-beta_r*Sr[i]/(Tr[i]) *F[i]*(1-s[i])-d_r*Sr[i]-rho*Sr[i])
        Ir[i+1]=Ir[i]+h*(beta_r*Sr[i]/Tr[i]*F[i]*(1-s[i])-(gamma_r+d_r)*Ir[i])
        Rr[i+1]=Rr[i]+h*(r_r*Rr[i]*(p-Tr[i]/K_r)+gamma_r*g_r*Ir[i]-d_r*Rr[i])
        Dr[i+1]=Dr[i]+h*(+gamma_r*(1-g_r)*Ir[i]+d_r*Tr[i])
        Vr[i+1]=Vr[i]+h*(Sr[i]*rho-d_r*Vr[i])
    
        
        Nf[i+1]=Nf[i]+h*(r_f*Nf[i]*(1-Nf[i]/K_f)+mu_f/Tr[i]*F[i]*(1-s[i]))
        F[i+1]=F[i]+h*((d_r+gamma_r*(1-g_r))*Ir[i]*Nf[i]-d_f*F[i])
        
        Sh[i+1]=Sh[i]+h*(-beta_h*Sh[i]*F[i]*s[i])
        Ih[i+1]=Ih[i]+h*(beta_h*Sh[i]*F[i]*s[i]-gamma_h*Ih[i])
        Rh[i+1]=Rh[i]+h*(gamma_h*Ih[i]*(1-mu_h))
        Dh[i+1]=Dh[i]+h*(gamma_h*Ih[i]*mu_h)
        
        
    Tr[N-1]=Sr[N-1]+Ir[N-1]+Rr[N-1]+Vr[N-1]
    
    
    #plt.fill_between(t,Dh+Ih+Rh,Dh+Sh+Ih+Rh,label='Susceptible')
    #plt.fill_between(t,Dh+Rh,Dh+Ih+Rh,label='Infected')
    #plt.fill_between(t,Dh,Dh+Rh,label='Recovered')
    #plt.fill_between(t,Dh,label='Dead')
    #plt.legend()
    #plt.show()
    
    
    #plt.plot(t*365,Dh,label='Modelled Deaths')
    #plt.plot(tlst,totaldeaths,label='Actual deaths')
    #plt.legend()
    #plt.show()
    
    
    #plt.plot(t,Sr,label='Sr')
    #plt.plot(t,Ir,label='Ir')
    #plt.plot(t,Rr,label='Rr')
    #plt.plot(t,Tr,label="Tr")

    #plt.plot(t,Vr,label="Vr")
    #plt.plot(t,Dr,label='D')
    #plt.legend()
    #plt.show()
    return(Dh[N-1])

loops=100
rho=np.linspace(0,0.2*365,loops)
D=np.zeros(loops)

for i in range(loops):
    D[i]=loop(rho[i])
    
plt.plot(rho/365,D)
plt.ylabel("Deaths")
plt.ylim(bottom=0)

plt.xlabel("Vaccination rate")
plt.title("Vaccinating rats in flea-based model")
plt.show()


    
  