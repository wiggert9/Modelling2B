# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:55:10 2022

@author: Wiggert
"""




import numpy as np
import matplotlib.pyplot as plt

S=0.99
I=1-S
R=0
h=0.01
T=200
N=int(T/h)
Sarr=np.zeros(N)
Iarr=np.zeros(N)
Rarr=np.zeros(N)
tarr=np.zeros(N)

R0=1.6666
gamma=0.2
beta=R0*gamma
t=0

def function(array):
    S,I,R = array[0],array[1],array[2]
    term1 = -beta*S*I
    term2 = beta*S*I-gamma*I
    term3 = gamma*I
    
    return np.array([term1,term2,term3])

def rungekatta(array):
    k1 = h*function(array)
    k2=h*function(array+0.5*k1)
    k3=h*function(array+0.5*k2)
    k4 = h*function(array+k3)
    
    return array+1/6 * (k1+2*k2+2*k3+k4)
    

for i in range(N):
    Sarr[i]=S
    Iarr[i]=I
    Rarr[i]=R
    tarr[i]=t
    array=np.array([S,I,R])
    newarray=rungekatta(array)
    S,I,R = newarray[0],newarray[1],newarray[2]
    t=t+h
    
    
    
        
    
plt.plot(tarr,Sarr,label='S')
plt.plot(tarr,Iarr,label='I')
plt.plot(tarr,Rarr,label='R')
plt.legend()
plt.show()


plt.fill_between(tarr,Iarr+Rarr,Sarr+Iarr+Rarr,label='Susceptible')
plt.fill_between(tarr,Rarr,Iarr+Rarr,label='Infected')
plt.fill_between(tarr,Rarr,label='Recovered')
plt.legend()
plt.show()



#starting april 25th


