# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:14:34 2022

@author: María Pérez
"""

import matplotlib.pyplot as plt

r=[]

C1=[]

C2=[]

with open("CorrelationFunction0.25.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        r.append(float(vals[0]))
        
        C1.append(float(vals[1]))
        
with open("CorrelationFunction0.50.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        C2.append(float(vals[1]))

plt.plot(r, C1, label='T/J=0.25')
plt.plot(r, C2, label='T/J=0.50')
plt.ylim([0,1.1])
plt.legend()
plt.grid(True)
plt.xlabel(r'$r$')
plt.ylabel(r'$C(r)$')
plt.title(r'Real part of the correlation function vs r. for two values of the temperature')
#plt.savefig('Correlation function.pdf')
