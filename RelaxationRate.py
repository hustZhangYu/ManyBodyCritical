# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 11:08:41 2023

perturbation calculation for  the relaxation rate

@author: YuZhang
"""

import numpy as np
import spin as sp
import random as rd


def GetMatrixSingleMu(Ev, Lk):
# =============================================================================
#     Ev :Eigenvectors 
#     Lk: jump operator
#     
# =============================================================================

    L=np.size(Ev,0)
    M=np.zeros([L,L])
    for m in range(L):
        PsiM=Ev[:,m]
        for n in range(L):
            PsiN=Ev[:,n]
            if m==n:
                Ma=np.dot(np.conj(PsiM.T),Lk)
                Mb=np.dot(Ma,Lk)
                Mc=np.dot(Mb,PsiN)
                M1=np.dot(np.conj(PsiM.T),Lk)
                M1=np.dot(M1,PsiN)
                M[m,n]=Mc-M1*np.conj(M1)
            else:
                M1=np.dot(np.conj(PsiM.T),Lk)
                M1=np.dot(M1,PsiN)
                M[m,n]=-M1*np.conj(M1)
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