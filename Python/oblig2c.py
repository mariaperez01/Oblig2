# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:06:13 2022

@author: María Pérez Martínez
"""
import matplotlib.pyplot as plt

T=[]

m16=[]

        
with open("m_x16forT0to1000100pointscambio.txt") as infile:
    infile.readline()
    
    lines= infile.readlines()
    
    for line in lines:
        
        vals = line.split()
        
        T.append(float(vals[0]))
        
        m16.append(float(vals[1]))


print("The value of <m> at T/J=0 is:", m16[0])

print("The value of <m> at T/J= inf is:", m16[-1])

plt.plot(T, m16)
plt.legend()
plt.grid(True)
plt.xlabel(r'$T/J$')
plt.ylabel(r'<$m$/Am$^-1$>')
plt.title(r'Average magnetization per site <$m$> vs. $T/J$ for $L$=16')
