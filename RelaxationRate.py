# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 11:08:41 2023

perturbation calculation for  the relaxation rate

@author: YuZhang
"""

import numpy as np
import spin as sp
import random as rd
import time 


def GetMatrixSingleMu1(Ev, Lk):
# =============================================================================
#     Ev :Eigenvectors 
#     Lk: jump operator
# =============================================================================
    print('ok1')
    L=np.size(Ev,0)
    M=np.zeros([L,L])
    for m in range(L):
        PsiM=Ev[:,m]
        PsiM1=np.conj(PsiM.T);
        print('m_begin')
        time1=time.time()
        for n in range(L):
            PsiN=Ev[:,n]
# =============================================================================
#             print(m,n)
# =============================================================================
            if m==n:
                Ma=np.dot(PsiM1,Lk)
                Mb=np.dot(Ma,Lk)
                Mc=np.dot(Mb,PsiN)
                M1=np.dot(PsiM1,Lk)
                M1=np.dot(M1,PsiN)
                M[m,n]=Mc-M1*np.conj(M1)
            else:
                M1=np.dot(PsiM1,Lk)
                M1=np.dot(M1,PsiN)
                M[m,n]=-M1*np.conj(M1)
        time2=time.time()
        print(time2-time1)
    print('ok2')
    return M


# =============================================================================
# def GetMatrixSingleMu(Ev, Lk):
# # =============================================================================
# #     Ev :Eigenvectors 
# #     Lk: jump operator
# # =============================================================================
#     timea=time.time()    
#     print('timea=',timea)    
#     L=np.size(Ev,0)
#     M=np.zeros([L,L])
#     
#     Ev1=np.conj(Ev.T);
#     Fterm0=np.diag(Ev1@Lk@Lk@Ev)
#     Fterm=np.diag(Fterm0)
#     
#     Sterm0=Ev1@Lk@Ev
#     Sterm=Sterm0*np.conj(Sterm0)
#     
#     M=Fterm-Sterm
#     timeb=time.time()  
#     print('timeb=',timeb-timea)   
#     return M
# =============================================================================


def GetMatrixSingleMu2(Ev, Lk):
 # =============================================================================
 #     Ev :Eigenvectors 
 #     Lk: jump operator
 # =============================================================================
     timea=time.time()    
     print('timea=',timea)    
     L=np.size(Ev,0)
     M=np.zeros([L,L])
     
     Ev1=np.conj(Ev.T);
     Fterm0=np.diag(Ev1@Lk@Lk@Ev)
     Fterm=np.diag(Fterm0)
     
     Sterm0=Ev1@Lk@Ev
     Sterm=Sterm0*np.conj(Sterm0)
     
     M=Fterm-Sterm
     timeb=time.time()  
     print('timeb=',timeb-timea)   
     return M


def GetMatrixSingleMu(Ev, Lk):
# =============================================================================
#     Ev :Eigenvectors 
#     Lk: jump operator
# =============================================================================
    timea=time.time()    
    print('timea=',timea)    
    L=np.size(Ev,0)
    
    print('a')
    M=np.zeros([L,L])
    
    print('b')
    Ev1=np.conj(Ev.T);
    
    print('c')
    M=np.diag(np.diag(Ev1@Lk@Lk@Ev))-(Ev1@Lk@Ev)*np.conj(Ev1@Lk@Ev)
    
    print('d')
    timeb=time.time()  
    print('timeb=',timeb-timea)   
    return M


def RelaxationRate(M):    
# =============================================================================
#     M: M perturbation matrix
# =============================================================================
    E,Ev=np.linalg.eig(M)
    E1=np.sort(np.real(E))
    return E1


if __name__=="__main__":
    x=1
    print(x)