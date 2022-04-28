# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:18:46 2022

@author: maari
"""
import matplotlib.pyplot as plt

T=[]

gamma8=[]

gamma16=[]

gamma32=[]

with open("m_x8forT0toT21000points.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        T.append(float(vals[0]))
        
        gamma8.append(float(vals[4])/8)
        
with open("m_x16forT0toT21000points.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        gamma16.append(float(vals[4])/16)

with open("m_x32forT0toT21000points.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        gamma32.append(float(vals[4])/32)

plt.plot(T, gamma8, label='L=8')
plt.plot(T, gamma16, label='L=16')
plt.plot(T, gamma32, label='L=32')
plt.legend()
plt.grid(True)
plt.xlim([0,2])
plt.xlabel(r'$T/J$')
plt.ylabel(r'$\Gamma$/L')
plt.title(r'Curves of $\Gamma$ vs. $T/J$ for differente values of $L$')
#plt.savefig('Correlation function.pdf')