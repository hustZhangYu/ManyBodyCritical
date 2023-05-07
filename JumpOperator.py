# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:03:50 2023

We define the Lindblad operator in spin systems

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

def SigmaZ_i(L,n,m):
     
    # m is the m_th site, no need to -1 (F**k python)
    # return the O operator
    # spin operator
    sigmax=sp.sigma_x()
    sigmay=sp.sigma_y()
    sigmaz=sp.sigma_z()
    I=sp.I()
    

    h2=np.zeros(L)
    h2[m-1]=1
    #Hamiltonian
    O=np.zeros([2**L,2**L])
    
    # disorder term
    for i in range(L):
        O=O+np.kron(np.eye(2**i),np.kron(h2[i]*sigmaz,np.eye(2**(L-i-1))))   
    
# =============================================================================
#     # fixed particle number
#     list1=nlabel(L,n)
#     La=len(list1)
#     O1=np.zeros([La,La])
#     for i in range(La):
#         for j in range(La):
#             O1[i,j]=O[list1[i],list1[j]]
# =============================================================================

    return O

def SigmaX_i(L,n,m):
     
    # m is the m_th site, no need to -1 (F**k python)
    # return the O operator
    # spin operator
    sigmax=sp.sigma_x()
    I=sp.I()
    

    h2=np.zeros(L)
    h2[m-1]=1

    #Hamiltonian
    O=np.zeros([2**L,2**L])
    
    # disorder term
    for i in range(L):
        O=O+np.kron(np.eye(2**i),np.kron(h2[i]*sigmax,np.eye(2**(L-i-1))))   
    

    return O

def SigmaY_i(L,n,m):
     
    # m is the m_th site, no need to -1 (F**k python)
    # return the O operator
    # spin operator
    sigmay=sp.sigma_y()
    I=sp.I()
    

    h2=np.zeros(L)
    h2[m-1]=1

    #Hamiltonian
    O=np.zeros([2**L,2**L])
    
    # disorder term
    for i in range(L):
        O=O+np.kron(np.eye(2**i),np.kron(h2[i]*sigmay,np.eye(2**(L-i-1))))   

    return O

    

if __name__=="__main__":
    X=1
    print(X)


