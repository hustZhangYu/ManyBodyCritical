# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:51:57 2023

Hamiltonian for the MBC region

@author: Yu Zhang
"""

import numpy as np
import spin as sp
import random as rd

def nlabel(L,n):
    """
    Parameters
    ----------
    L : int
        size of the system.
    n : int
        numbers of particle
    Returns
    -------
    sequence of the basis with the particle number
    """
    No=[[1,0],[0,0]]
    N=np.zeros([2**L,2**L])
    for i in range(L):
        N=N+np.kron(np.eye(2**i),np.kron(No,np.eye(2**(L-i-1))))  
    list1=[]
    for i in range(2**L):
        if N[i,i]==n:
            list1.append(i)
    return list1



def H1(L,mu,delta,V,U,n):
    
    # spin operator
    sigmax=sp.sigma_x()
    sigmay=sp.sigma_y()
    sigmaz=sp.sigma_z()
    I=sp.I()
    
    # constant
    omega=(np.sqrt(5)-1)/2
    
    # hopping term/ disorder/ interaction cofficents 
    # hopping term
    h1=(np.ones(L-1)+[mu*np.cos(2*np.pi*omega*(j+1/2)+delta) for j in range(L-1)])/2
    # disorder
    h2=[V*np.cos(2*np.pi*omega*j+delta) for j in range(L)]
    # interaction
    h3=[U for j in range(L-1)]
    
    
    #Hamiltonian
    H0=np.zeros([2**L,2**L])
    
    # hopping term and interaction
    for i in range(L-1):
        H0=H0+np.kron(np.eye(2**i),np.kron(h1[i]*(np.kron(sigmax,sigmax)+np.kron(sigmay,sigmay))+h3[i]*np.kron((sigmaz+I)/2,(sigmaz+I)/2),np.eye(2**(L-i-2))))
    # disorder term
    for i in range(L):
        H0=H0+np.kron(np.eye(2**i),np.kron(h2[i]*(sigmaz+I)/2,np.eye(2**(L-i-1))))   
    
    # fixed particle number
    list1=nlabel(L,n)
    La=len(list1)
    H=np.zeros([La,La])
    for i in range(La):
        for j in range(La):
            H[i,j]=H0[list1[i],list1[j]]

    return H




if __name__=="__main__":
    X=1
    print(X)
