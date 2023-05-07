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
import time


def MbcR(L,n,mu,delta,V,U):
    
# =============================================================================
#     
#     L:system size 
#     n: particle number 
#     delta: random phase 
#     V:the strength of the disorder
#     U: nearest neighbour interaction strength 
#     E1: all the perturbation energy 
# ========================================= ====================================
    
    m=1
    H=MbcH.H1(L,mu,delta,V,U,n)
    La=JP.SigmaZ_i(L,n,m)
    Lb=JP.SigmaX_i(L,n,m)
    Lc=JP.SigmaY_i(L,n,m)
    time1=time.time()
    print("time1=" ,time1)
    E,Ev=np.linalg.eig(H)
    M=RR.GetMatrixSingleMu(Ev, La)+RR.GetMatrixSingleMu(Ev, Lb)+RR.GetMatrixSingleMu(Ev, Lc)
    M1=RR.GetMatrixSingleMu2(Ev, La)+RR.GetMatrixSingleMu2(Ev, Lb)+RR.GetMatrixSingleMu2(Ev, Lc)
    print(M-M1)
    time2=time.time()
    print(time2-time1)
# =============================================================================
#     M=RR.GetMatrixSingleMu(Ev, La)
# =============================================================================
# =============================================================================
#     M1=RR.GetMatrixSingleMu1(Ev, La)+RR.GetMatrixSingleMu1(Ev, Lb)+RR.GetMatrixSingleMu1(Ev, Lc)
# =============================================================================

    E1=RR.RelaxationRate(M)
    E1=E1*np.size(H,0)
    return E1


def ChangeW(x):  
    
# =============================================================================
#     W: from 10^(-1) to 10^2
#     ls: length of the data
# =============================================================================
    
    L=12; n=int(L/2); mu=1.5; delta=2*np.pi*rd.random(); U=1;
    
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

def ChangeWTest():  
    
# =============================================================================
#     W: from 10^(-1) to 10^2
#     ls: length of the data
# =============================================================================
    
    L=12; n=int(L/2); mu=1.5; delta=2*np.pi*rd.random(); U=1;
        
    ls=41;
    W_all=[10**(-1+i*0.1) for i in range(ls)]
    W_all=[0.5*i for i in range(ls)]
    data=np.zeros(ls)
    print(delta)    
    for i in range(ls):
        W=W_all[i]
        t1=time.time()
        E1=MbcR(L,n,mu,delta,W,U)
# =============================================================================
#         E1=(1,2)          
# =============================================================================
        data[i]=E1[1]
        t2=time.time()
        print(t2-t1)
    plt.plot(W_all,np.log(data))
    return data

def ChangeW10(x):  
    
# =============================================================================
#     W: from 10^(-1) to 10^2
#     ls: length of the data
# =============================================================================
    
    L=10; n=int(L/2); mu=1.5; delta=2*np.pi*rd.random(); U=1;
    
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

def ChangeW12(x):  
    
# =============================================================================
#     W: from 10^(-1) to 10^2
#     ls: length of the data
# =============================================================================
    
    L=12; n=int(L/2); mu=1.5; delta=2*np.pi*rd.random(); U=1;
    
    ls=30;
    W_all=[10**(-1+i*0.1) for i in range(ls)]
    data=np.zeros(ls)
    
    for i in range(ls):
        W=W_all[i]
        time1=time.time()
        print(time1)
        E1=MbcR(L,n,mu,delta,W,U)
# =============================================================================
#         E1=(1,2)          
# =============================================================================
        data[i]=E1[1]
        time2=time.time()
        print(time2-time1)
    return data

def ChangeW14(x):  
    
# =============================================================================
#     W: from 10^(-1) to 10^2
#     ls: length of the data
# =============================================================================
    
    L=10; n=int(L/2); mu=1.5; delta=2*np.pi*rd.random(); U=1;
    
    ls=30;
    W_all=[10**(-1+i*0.1) for i in range(ls)]
    data=np.zeros(ls)
    
    for i in range(ls):
        W=W_all[i]
        time1=time.time()
        print(time1)
        E1=MbcR(L,n,mu,delta,W,U)
# =============================================================================
#         E1=(1,2)          
# =============================================================================
        data[i]=E1[1]
        time2=time.time()
        print(time2-time1)
    return data

def RealxationRatioSample():
    sample=1000
    L=10
    pool=multiprocessing.Pool(processes = 52)
    res=pool.map(ChangeW10, range(sample))
    pool.close()    
    pool.join()
    df=pd.DataFrame(res)
    df.to_csv('L{}.csv'.format(L),index=False,header=False)
    sample=1000
    L=12
    pool=multiprocessing.Pool(processes = 52)
    res=pool.map(ChangeW12, range(sample))
    pool.close()    
    pool.join()
    df=pd.DataFrame(res)
    df.to_csv('L{}.csv'.format(L),index=False,header=False)
    sample=300
    L=14
    pool=multiprocessing.Pool(processes = 52)
    res=pool.map(ChangeW14, range(sample))
    pool.close()    
    pool.join()
    df=pd.DataFrame(res)
    df.to_csv('L{}.csv'.format(L),index=False,header=False)


def RealxationRatioSample1():
    
    ChangeW14(1)

    
# =============================================================================
#     sample=10
#     time1=time.time()
#     print(time1)
#     L=12
#     pool=multiprocessing.Pool(processes = 3)
#     res=pool.map(ChangeW12, range(sample))
#     pool.close()    
#     pool.join()    
#     time2=time.time()
#     print(time2)
#     print(time2-time1)
# =============================================================================


def RaSampleDistribution(Data):    
# =============================================================================
#     We plot the distribution of Relaxation rate after taking some samples.
#     Input： Data (Relaxation ratio from different samples)
# =============================================================================
    
    x=1

    return x


if __name__=="__main__":
    RealxationRatioSample1()
    
# =============================================================================
#     L=6; n=3; mu=1.5; delta=2*np.pi*rd.random(); U=1;
#     V=10;
#     m=1
#     H=MbcH.H1(L,mu,delta,V,U,n)
       #     La=JP.SigmaZ_i(L,n,m)
#     Lb=JP.SigmaX_i(L,n,m)
#     Lc=JP.SigmaY_i(L,n,m)
#     E,Ev=np.linalg.eig(H)
#     M=RR.GetMatrixSingleMu(Ev, La)+RR.GetMatrixSingleMu(Ev, Lb)+RR.GetMatrixSingleMu(Ev, Lc)
#     M1=RR.GetMatrixSingleMu1(Ev, La)+RR.GetMatrixSingleMu1(Ev, Lb)+RR.GetMatrixSingleMu1(Ev, Lc)
#     print(M)
#     print(M1)
# =============================================================================
