#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of utility functions for different parts of the code

Created on Fri Dec  6 16:08:35 2024
@author: mmeb1
"""
import numpy as np
from scipy.interpolate import interp1d

#%% Functions for different iceberg area profiles wrt depth
def iceberg_fun_exponential1(nu,hgl,z): # we take all 3 arguments for consistency when assigning function to an object attribute
    return np.exp(-z/100)

def iceberg_fun_exponential2(nu,hgl,z):
    return (nu/hgl)*np.exp(nu*-z/hgl)/(1-np.exp(-nu))

#%% Functions for fjord hypsometry
def hypsometry_idealised(width,z): #TODO: implement idealised hypsometry
    pass

#%% Other miscellaneous functions
# quick function to find nearest point, especially useful for the time axis
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))



#%% binning forcings to a fjord's vertical profile

def bin_forcings(f, H, t):
    """
    BIN_FORCINGS Puts forcings on model layers and time steps.
    [Ts, Ss, Qsg] = BIN_FORCINGS(f, H, t) calculates the mean value of 
    shelf temperature and salinity over the depths of each model layer. It 
    then interpolates both the resulting shelf profiles and subglacial 
    discharge onto the model time steps.
    
    Parameters:
    - f: Dictionary containing forcing data, with keys 'Ts', 'Ss', 'zs', 'ts', 'Qsg', 'tsg'.
    - H: 1D array of layer thicknesses for the model.
    - t: 1D array of model time steps.
    
    Returns:
    - Ts: 2D array of temperature profiles for each model layer and time step.
    - Ss: 2D array of salinity profiles for each model layer and time step.
    - Qsg: 2D array of subglacial discharge for each plume and time step.
    """
    # Replace NaN entries in the shelf profiles by the last valid value
    f['Ts'] = f['Ts'].ffill(dim='depth')
    f['Ss'] = f['Ss'].ffill(dim='depth')
    #f['zs'] = f['Ts'].depth

    # First, put shelf forcing on model layers
    # Add the layer boundaries H0 into the vector of shelf z-values
    z0 = np.unique(np.sort(np.concatenate(([0], f['zs'], np.cumsum(H)))) )

    # Interpolate the shelf temperature and salinity profiles onto the new grid z0    
    T0 = np.zeros((len(z0), len(f['ts'])))  # Assuming f['Ts'] has multiple time steps
    S0 = np.zeros((len(z0), len(f['ts'])))
    for k in range(len(f['ts'])):
        T0[:,k] = np.interp(z0,f['zs'],f['Ts'][:,k])
        S0[:,k] = np.interp(z0,f['zs'],f['Ss'][:,k])


    # Calculate shelf temperature and salinity (Ts0 and Ss0) on model layers
    ints = np.concatenate(([0], np.cumsum(H)))
    Sz = np.zeros((len(H), len(f['ts'])))  # Assuming f['Ts'] has multiple time steps
    Tz = np.zeros((len(H), len(f['ts'])))
    
    for k in range(len(ints) - 2):
        # Find the boundaries of the layer in the shelf grid z0
        inds = np.where((z0 >= ints[k]) & (z0 <= ints[k + 1]))[0]
        # Average the temperature and salinity profiles over this layer
        Sz[k, :] = np.trapz(S0[inds, :], x=z0[inds],axis=0) / H[k]
        Tz[k, :] = np.trapz(T0[inds, :], x=z0[inds],axis=0) / H[k]

    # Second, put forcings on model time steps
    # Interpolate shelf conditions to model time steps
    Ss = np.zeros((len(H), len(t)))
    Ts = np.zeros((len(H), len(t)))
    for k in range(len(H)):
        Ss[k,:] = np.interp(t, f['ts'], Sz[k,:])
        Ts[k,:] = np.interp(t, f['ts'], Tz[k,:])
    

    # Subglacial discharge
    Qsg = np.zeros((len(f['Qsg']), len(t)))
    for j in range(f['Qsg'].shape[0]):
        Qsg[j, :] = interp1d(f['tsg'], f['Qsg'][j, :], kind='linear', fill_value="extrapolate")(t)
    
    return Ts, Ss, Qsg