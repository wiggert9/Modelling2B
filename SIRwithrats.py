# -*- coding: utf-8 -*-
"""
Created on Wed May 25 13:04:19 2022

@author: Wiggert
"""
import numpy as np
import matplotlib.pyplot as plt



#starting april 25th
deathsperweek=[2,0,9,3,14,17,43,112,168,267,470,725,1089,1843,2010,2817,3880,4237,6102,6988,6544,7165,5533,4929,4327,2665,1421,1031,1414,1050,652,333,210,243,281]
totaldeaths=np.cumsum(deathsperweek)
#end december 19th

#plt.plot(deathsperweek)
tlst=np.linspace(0,34,35)*7


plt.show()

#plt.plot(tlst,totaldeaths)
#plt.xlabel("Time in days(since April 25th)")
#plt.ylabel("Total deaths")
plt.show()

population = 460000

Ih=0/population
Sh=1-Ih
Rh=0
Dh=0

Ir=1/population
Sr=1-Ir
Rr=0
Dr=0

h=0.01
T=238
N=int(T/h)
Sharr=np.zeros(N)
Iharr=np.zeros(N)
Rharr=np.zeros(N)
Dharr=np.zeros(N)

Srarr=np.zeros(N)
Irarr=np.zeros(N)
Rrarr=np.zeros(N)
Drarr=np.zeros(N)
tarr=np.zeros(N)

R0=3.2
gamma_1=1/26
beta_1=R0*gamma_1
gamma_2=1/5.15
beta_2=R0*gamma_2


beta=R0*gamma_1
t=0

mu_1=0.33
mu_2=0.9


for i in range(N):
    Sharr[i]=Sh
    Iharr[i]=Ih
    Rharr[i]=Rh
    Dharr[i]=Dh
    tarr[i]=t
    Shn,Ihn,Rhn,Dhn=Sh,Ih,Rh,Dh
    
    Srarr[i]=Sr
    Irarr[i]=Ir
    Rrarr[i]=Rr
    Drarr[i]=Dr
    tarr[i]=t
    Srn,Irn,Rrn,Drn=Sr,Ir,Rr,Dr
    
    
    Sh=Shn+h*(-beta_1*Shn*Irn)
    Ih=Ihn+h*(beta_1*Shn*Irn -gamma_1*Ihn)
    Rh=Rhn+h*((1-mu_1)*gamma_1*Ihn)
    Dh=Dhn+h*((mu_1*gamma_1*Ihn))
    t=t+h
    
    Sr=Srn+h*(-beta_2*Srn*Irn)
    Ir=Irn+h*(beta_2*Srn*Irn -gamma_2*Irn)
    Rr=Rrn+h*((1-mu_2)*gamma_2*Irn)
    Dr=Drn+h*((mu_2*gamma_2*Irn))
    
        

plt.plot(tarr,Sharr,label='Susceptible')
plt.plot(tarr,Iharr,label='Infected')
plt.plot(tarr,Rharr,label='Recovered')
plt.plot(tarr,Dharr,label='Dead')
plt.legend()
plt.xlabel("Days (since April 25th)")
plt.ylabel("Proportion in each compartment")
plt.title("Time evolution of distribution of humans over the SIRD comparments")
plt.show()


plt.fill_between(tarr,Dharr+Iharr+Rharr,Dharr+Sharr+Iharr+Rharr,label='Susceptible')
plt.fill_between(tarr,Dharr+Rharr,Dharr+Iharr+Rharr,label='Infected')
plt.fill_between(tarr,Dharr,Dharr+Rharr,label='Recovered')
plt.fill_between(tarr,Dharr,label='Dead')
plt.legend()
plt.ylabel("Proportion in each compartment")
plt.xlabel("Days (since April 25th)")
plt.title("SIRD model with rats with total population of humans shown divided over each compartment")
plt.show()


plt.plot(tarr,Dharr*population,label='Model with rats')
plt.plot(tlst,totaldeaths,label='Actual deaths')
plt.legend()
plt.xlabel("Days (since April 25th)")
plt.title("Total deaths predicted by SIRD model with rats compared to actual deaths")
plt.ylabel("Number of deaths")
plt.show()


plt.plot(tarr,Srarr,label='Susceptible')
plt.plot(tarr,Irarr,label='Infected')
plt.plot(tarr,Rrarr,label='Recovered')
plt.plot(tarr,Drarr,label='Dead')
plt.legend()
plt.xlabel("Days (since April 25th)")
plt.title("Time evolution of distribution of rats over the SIRD comparments")
plt.show()


plt.fill_between(tarr,Drarr+Irarr+Rrarr,Drarr+Srarr+Irarr+Rrarr,label='Susceptible')
plt.fill_between(tarr,Drarr+Rrarr,Drarr+Irarr+Rrarr,label='Infected')
plt.fill_between(tarr,Drarr,Drarr+Rrarr,label='Recovered')
plt.fill_between(tarr,Drarr,label='Dead')
plt.legend()
plt.ylabel("Proportion in each compartment")
plt.xlabel("Days (since April 25th)")
plt.title("SIRD model with rats with total population of rats shown divided over each compartment")
plt.show()
