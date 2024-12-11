#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of utility functions for different parts of the code

Created on Fri Dec  6 16:08:35 2024
@author: mmeb1
"""
import numpy as np

#%% Functions for different iceberg area profiles wrt depth
def iceberg_fun_exponential1(nu,hgl,z): # we take all 3 arguments for consistency when assigning function to an object attribute
    return np.exp(-z/100)

def iceberg_fun_exponential2(nu,hgl,z):
    return (nu/hgl)*np.exp(nu*z/hgl)/(1-np.exp(-nu))

#%% Functions for fjord hypsometry
def hypsometry_idealised(width,z): #TODO: implement idealised hypsometry
    pass

#%% Other miscellaneous functions
# quick function to find nearest point, especially useful for the time axis
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))
