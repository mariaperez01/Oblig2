# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:14:34 2022

@author: María Pérez
"""

import matplotlib.pyplot as plt
import numpy as np

T1=0.25

T2 =0.5

r=[]

N=16

C1=[]

C2=[]

lambda0=2+np.exp(1/T1)

lambda2=np.exp(1/T1)-1

lambda0T2=2+np.exp(1/T2)

lambda2T2=np.exp(1/T2)-1

with open("CorrelationFunction0.25.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        r.append(float(vals[0]))
        

for i in range(0,len(r)):
    
    C1.append(((lambda0**(r[i]))*lambda2**(N-r[i])+(lambda2**(N))+(lambda2**(r[i]))*(lambda0**(N-r[i])))/(lambda0**N+2*lambda2**N))
    
    C2.append(((lambda0T2**(r[i]))*lambda2T2**(N-r[i])+(lambda2T2**(N))+(lambda2T2**(r[i]))*(lambda0T2**(N-r[i])))/(lambda0T2**N+2*lambda2T2**N))


plt.plot(r, C1, label='T/J=0.25')
plt.plot(r, C2, label='T/J=0.50')
plt.legend()
plt.grid(True)
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
plt.title(r'Real part of the correlation function vs r. for two values of the temperature.')
#plt.savefig('Correlation function.pdf')
