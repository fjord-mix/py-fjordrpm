#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 17:45:29 2024

@author: mmeb1
"""
import numpy  as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import datetime as dt

#%% time vector for forcings
#time_begin = pd.to_datetime('2010-01-01')
#time_end   = pd.to_datetime('2013-01-01')
time_begin = dt.datetime(2010,1,1)
time_end   = dt.datetime(2013,1,1)

dt_model   = 0.2
ts         = np.array([time_begin,time_end])

#%% shelf
# depth vector for shelf forcing
zs = np.linspace(0,800,40)

# create shelf T and S
t_shelf = 3*np.ones([len(zs),len(ts)])
s_shelf = 34*np.ones([len(zs),len(ts)])

#%% Subglacial discharge forcing
tsg = np.arange(time_begin.to_julian_date(),time_end.to_julian_date(),dt_model) # create time axis

# here use idealised seasonal gaussian peaked at julian day 200
n_plumes = 1 # choose number of plumes here
i_plumes = [i for i in range(n_plumes)] 

# set up GL depths
gl_depth = 800.0
gl_depth_all_plumes = np.ones(n_plumes)*gl_depth



qsg = 300*np.exp(-((tsg % 365-200)/30)**2); # subglacial discharge on tsg
qsg_all_plumes = np.repeat(qsg[np.newaxis,:],n_plumes,axis=0)

#%% set up output datasets
ds_shelf = xr.Dataset(coords=dict(depth=(['depth'],zs),time=(['time'],ts)), \
                      attrs=dict(run_name="example1_subglacial_discharge"))
    
ds_glacier = xr.Dataset(coords=dict(plume=(['plume'],i_plumes),time=(['time'],tsg)), \
                      attrs=dict(run_name="example1_subglacial_discharge"))
    
    
ds_Ts  = xr.DataArray(t_shelf,coords={'depth':zs,'time':ts})
ds_Ss  = xr.DataArray(s_shelf,coords={'depth':zs,'time':ts})
ds_Qsg = xr.DataArray(qsg_all_plumes,coords={'plume':i_plumes,'time':tsg})
ds_Hgl = xr.DataArray(gl_depth_all_plumes,coords={'plume':i_plumes})

ds_shelf['thetao'] = ds_Ts
ds_shelf['so'] = ds_Ss
ds_glacier['Qsg'] = ds_Qsg
ds_glacier['Hgl'] = ds_Hgl

#%% saving...
ds_shelf.to_netcdf('example1_shelf_forcings.nc')
ds_glacier.to_netcdf('example1_glacier_forcings_np'+str(n_plumes)+'.nc')
