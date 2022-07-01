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
T=100
N=int(T/h)

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
T_r=np.zeros(N)

s=np.zeros(N)
t=np.zeros(N)

Sh[0]=population
Sr[0]=55000
Rr[0]=30000
Ir[0]=500






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

for i in range (N-1):
    T_r[i]=Sr[i]+Ir[i]+Rr[i]
    s[i]=np.exp(-a*T_r[i])
    t[i+1]=t[i]+h

    Sr[i+1]=Sr[i]+h*(r_r*Sr[i]*(1-T_r[i]/K_r)+r_r*Rr[i]*(1-p)-beta_r*Sr[i]/(T_r[i]) *F[i]*(1-s[i])-d_r*Sr[i])
    Ir[i+1]=Ir[i]+h*(beta_r*Sr[i]/T_r[i]*F[i]*(1-s[i])-(gamma_r+d_r)*Ir[i])
    Rr[i+1]=Rr[i]+h*(r_r*Rr[i]*(p-T_r[i]/K_r)+gamma_r*g_r*Ir[i]-d_r*Rr[i])
    Dr[i+1]=Dr[i]+h*(+gamma_r*(1-g_r)*Ir[i]+d_r*T_r[i])
    
    Nf[i+1]=Nf[i]+h*(r_f*Nf[i]*(1-Nf[i]/K_f)+mu_f/T_r[i]*F[i]*(1-s[i]))
    F[i+1]=F[i]+h*((d_r+gamma_r*(1-g_r))*Ir[i]*Nf[i]-d_f*F[i])
    
    Sh[i+1]=Sh[i]+h*(-beta_h*Sh[i]*F[i]*s[i])
    Ih[i+1]=Ih[i]+h*(beta_h*Sh[i]*F[i]*s[i]-gamma_h*Ih[i])
    Rh[i+1]=Rh[i]+h*(gamma_h*Ih[i]*(1-mu_h))
    Dh[i+1]=Dh[i]+h*(gamma_h*Ih[i]*mu_h)
  
T_r[N-1]=Sr[N-1]+Ir[N-1]+Rr[N-1]


plt.plot(t,Sr,label='S')
plt.plot(t,Ir,label='I')
plt.plot(t,Rr,label='R')
#plt.plot(t,Dr,label='D')
plt.xlabel("time (years)")
plt.ylabel("number of rats")
plt.title("SIRD model for rats")
plt.legend()
plt.show()

plt.plot(t,Sh,label='S')
plt.plot(t,Ih,label='I')
plt.plot(t,Rh,label='R')
plt.plot(t,Dh,label='D')
plt.title("SIRD model for humans")
plt.legend()
plt.show()
'''


plt.plot(t,T_r*Nf, label='Total number of fleas')
plt.plot(t,F, label='Number of free infectious fleas')
plt.title("")
plt.legend()
plt.show()



plt.fill_between(t,Dh+Ih+Rh,Dh+Sh+Ih+Rh,label='Susceptible')
plt.fill_between(t,Dh+Rh,Dh+Ih+Rh,label='Infected')
plt.fill_between(t,Dh,Dh+Rh,label='Recovered')
plt.fill_between(t,Dh,label='Dead')
plt.title("Different way of visualizing SIRD model humans")
plt.legend()
plt.show()


plt.plot(t*365,Dh,label='Modelled Deaths')
plt.plot(tlst,totaldeaths,label='Actual deaths')
plt.title("Model compared to data")
plt.legend()
plt.show()




plt.fill_between(t,Dr+Ir+Rr,Dr+Sr+Ir+Rr,label='Susceptible')
plt.fill_between(t,Dr+Rr,Dr+Ir+Rr,label='Infected')
plt.fill_between(t,Dr,Dr+Rr,label='Recovered')
plt.fill_between(t,Dr,label='Dead')
plt.title("Different way of visualizing SIRD model rats")
plt.legend()
plt.show()
  
'''
  
    
  