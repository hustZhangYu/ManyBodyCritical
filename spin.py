# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:00:45 2022

@author: Lenovo
"""

import numpy as np




# spin operator
def sigma_x():
    sigmax=np.array([[0,1],[1,0]])
    return sigmax
def sigma_y():
    sigmay=np.array([[0,-1j],[1j,0]])
    return sigmay
def sigma_z():
    sigmaz=np.array([[1,0],[0,-1]])
    return sigmaz
def I():
    I=np.array([[1,0],[0,1]])
    return I
def sigma_plus():
    sigma_plus=(sigma_x()+1j*sigma_y())/2
    return sigma_plus
def sigma_minus():
    sigma_minus=(sigma_x()-1j*sigma_y())/2
    return sigma_minus





if __name__=="__main__":
    a=sigma_plus()
    print(a)
 