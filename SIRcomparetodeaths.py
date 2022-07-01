# -*- coding: utf-8 -*-
"""
Created on Wed May 25 12:11:18 2022

@author: Wiggert
"""



import numpy as np
import matplotlib.pyplot as plt



#starting april 25th
deathsperweek=[2,0,9,3,14,17,43,112,168,267,470,725,1089,1843,2010,2817,3880,4237,6102,6988,6544,7165,5533,4929,4327,2665,1421,1031,1414,1050,652,333,210,243,281]
totaldeaths=np.cumsum(deathsperweek)
#end december 19th

plt.plot(deathsperweek)
tlst=np.linspace(0,34,35)*7


plt.show()

plt.plot(tlst,totaldeaths)
plt.xlabel("Time in days(since April 25th)")
plt.ylabel("Total deaths")
plt.show()

population = 460000

I=10/population
S=1-I
R=0
D=0
h=0.01
T=300
N=int(T/h)
Sarr=np.zeros(N)
Iarr=np.zeros(N)
Rarr=np.zeros(N)
Darr=np.zeros(N)
tarr=np.zeros(N)

R0=3.2
gamma=1/(2.5+4.3) #Infectious period plus latent period 
beta=R0*gamma
t=0

mu=0.3

R0=3.2
gamma=1/(26)
beta=R0*gamma

Imax=I

for i in range(N):
    Sarr[i]=S
    Iarr[i]=I
    Rarr[i]=R
    Darr[i]=D
    tarr[i]=t
    Sn,In,Rn,Dn=S,I,R,D
    
    S=Sn+h*(-beta*Sn*In)
    I=In+h*(beta*Sn*In -gamma*In)
    R=Rn+h*((1-mu)*gamma*In)
    D=Dn+h*((mu*gamma*In))
    t=t+h
    if I>Imax:
        Imax=I
        Sturn=S
        k=i
    
    
    
        
    
plt.plot(tarr,Sarr,label='Susceptible')
plt.plot(tarr,Iarr,label='Infected')
plt.plot(tarr,Rarr,label='Recovered')
plt.plot(tarr,Darr,label='Dead')
plt.legend()
plt.title("SIRD model for bubonic plague")
plt.xlabel("Days (since April 25th)")
plt.ylabel("Proportion in each compartment")
plt.show()


plt.fill_between(tarr,Darr+Iarr+Rarr,Darr+Sarr+Iarr+Rarr,label='Susceptible')
plt.fill_between(tarr,Darr+Rarr,Darr+Iarr+Rarr,label='Infected')
plt.fill_between(tarr,Darr,Darr+Rarr,label='Recovered')
plt.fill_between(tarr,Darr,label='Dead')
plt.legend()
plt.xlabel("Days (since April 25th)")
plt.ylabel("Proportion in each compartment")
plt.title("SIRD model with total population shown divided over each compartment")
plt.show()


plt.plot(tarr,Darr*population,label='Model')
plt.plot(tlst,totaldeaths,label='Actual')
plt.legend()
plt.title("Total deaths predicted by model compared to actual deaths")
plt.xlabel("Days (since April 25th)")
plt.ylabel("Number of deaths")
plt.show()


