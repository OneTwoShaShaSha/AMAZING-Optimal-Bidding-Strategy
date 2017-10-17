#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 20:29:22 2017

@author: ayushbasu
"""

import numpy as np
from scipy.optimize import linprog

#############################################################
"""Program Constants"""
#N = Number of Prices
N = 50


#Delay Time
td = 100

#Necessary Working Time
W = 50


#Reboot Period
tr = 0.0833333


#Checkpointing Period
tc = 0.0833333


#SPOT PRICE ARRAY
E = np.zeros(N) 
for i in range(1,N+1):
    E[i-1] = i 
    
#This gives us a price array with a mean price of 50 and variation of 10

#For now we shall work with a random transition probability matrix
P = np.zeros(shape=(N,N))

for i in range(0,N):
    for j in range(0,N-1):
        P[i][j]=0.1
    P[i][N-1] = 0.51
     



#The above code ensures that it is a stochastic probability matrix by summing each row to 1. 
    

#############################################################
"""Defining the Basis Functions"""

#Computation Time Available in A State f
def f(i):
    if(i in range(1,N+1)):
        return(float(1))
    elif(i in range(N+1, 2*N+1)):
        return(float(1 - tr))
    elif(i in range(4*N+1, 5*N +1)):
        return(float(1-tc))
    elif(i in range(5*N+1, 6*N +1)):
        return(float(1-tc-tr))
    else:
        return(float(0))
    

#The Delta Kronecker Function
def delta(i,j):
    if ((i-j)%(4*N)==0):
        return 1.0
    else:
        return 0
    
##############################################################
"""Defining the Constraint Matrices"""

#Defining c or the objective functions coefficient
c = np.zeros(8*N)
for i in range(0,8*N):
    c[i]= E[(i%N)]*td
     
#print(c)

A = np.zeros(shape=(1,8*N))
#Defining A_ub or the inequality constraints  coefficient matrix
for i in range(0,8*N):
    A[0][i] = -f(i+1)
#Defining b_ub
b = [-W/td]

#print(A)
#print(b)


#Defining A_eq or the equality constraints coefficient matrix
A1 = np.zeros(shape=(8*N+1, 8*N))
b1 = np.zeros(8*N+1)
for i in range(0,8*N):
    A1[0][i] = 1
    b1[0]=1

for i in range(1, 8*N+1):
    for j in range(0,8*N):
        A1[i][j]= delta(i,j) - P[j%N][(i-1)%N]
    b1[i] = 0
          
#print(A1)
#print(b1)


#Bounds
x = [(0.0,None)]
x1 = [(0.0,None)]
for i in range(1,8*N):
    x = x+x1
x = tuple(x)

#print(x)


    
      
res = linprog(c,A_ub = A, b_ub = b, A_eq = A1, b_eq = b1, bounds = x, options ={"disp": True})

print(res)
       
    


    

    

