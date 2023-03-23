# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:28:20 2023

Main function for calculating the Relaxation rate of the MBC region.

@author: Yu Zhang 
"""


import numpy as np
import spin as sp
import random as rd
import RelaxationRate as RR
import JumpOperator as JP
import MbcH
import multiprocessing
from  matplotlib import pyplot as plt
import pandas as pd


def MbcR(L,n,mu,delta,V,U):
    
# =============================================================================
#     
#     L:system size 
#     n: particle number 
#     delta: random phase 
#     V:the strength of the disorder
#     U: nearest neighbour interaction strength 
#     E1: all the perturbation energy 
# =============================================================================
    
    m=1
    H=MbcH.H1(L,mu,delta,V,U,n)
    L=JP.SigmaZ_i(L,n,m)
    E,Ev=np.linalg.eig(H)
    M=RR.GetMatrixSingleMu(Ev, L)
    E1=RR.RelaxationRate(M)
    
    return E1


def ChangeW(x):  
    
# =============================================================================
#     W: from 10^(-1) to 10^2
#     ls: length of the data
# =============================================================================
    
    L=12; n=6; mu=1.5; delta=2*np.pi*rd.random(); U=1;
    
    ls=30;
    W_all=[10**(-1+i*0.1) for i in range(ls)]
    data=np.zeros(ls)
    
    for i in range(ls):
        W=W_all[i]
        E1=MbcR(L,n,mu,delta,W,U)
# =============================================================================
#         E1=(1,2)          
# =============================================================================
        data[i]=E1[1]
    return data

    
def ChangeWSample():   
# =============================================================================
#     multiprocess to take samples
#     In windows system, the multiprocess.Pool can't be written in function 
# =============================================================================
    sample=100
    pool=multiprocessing.Pool(processes = 3)
    res=pool.map(ChangeW, range(sample))
    pool.close()
    pool.join()
    df=pd.DataFrame(res)
    df.to_csv('L12.csv',index=False)



if __name__=="__main__":
    sample=100
    sequence=[i for i in range(sample)]
    pool=multiprocessing.Pool(processes = 3)
    res=pool.map(ChangeW, range(sample))
    pool.close()
    pool.join()
    df=pd.DataFrame(res)
    df.to_csv('L12.csv',index=False)
       