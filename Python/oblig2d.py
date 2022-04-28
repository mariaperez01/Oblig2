# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:14:34 2022

@author: maari
"""
import matplotlib.pyplot as plt

T=[]

m216=[]

        
with open("m_x16forT0to1000100pointscambio.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        T.append(float(vals[0]))
        
        m216.append(float(vals[2]))
        
print("The value of <|m|^2> at T/J=0 is:", m216[0])

print("The value of <|m|^2> at T/J= inf is:", m216[-1])

plt.plot(T, m216)
plt.legend()
plt.grid(True)
plt.xlabel(r'$T/J$')
plt.ylabel(r'<|$m$|$^2$>')
plt.title(r'Average magnetization squared per site <|$m$|$^2$> vs. $T/J$ for $L$=16')
