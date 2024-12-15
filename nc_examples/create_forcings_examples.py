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

example_number = 1

# depth vector for shelf forcing
zs = np.linspace(0,800,40)


# time vector for forcings
time_begin = pd.to_datetime('2010-01-01')
t_begin    = 0.0
dt_model   = 0.2

match example_number:
    case 1:
        # time axis
        t_end = 3*365
        ts    = np.array([t_begin,t_end]) # time vector for shelf forcing

        # shelf forcing
        t_shelf = 3*np.ones([len(zs),len(ts)])  # temperature
        s_shelf = 34*np.ones([len(zs),len(ts)]) # salinity
        
        # Subglacial discharge forcing
        tsg = np.arange(t_begin,t_end,dt_model)     # create time axis
        qsg = 300*np.exp(-((tsg % 365-200)/30)**2); # subglacial discharge on tsg (idealised seasonal gaussian peaked at julian day 200)
    case 2: 
        # time axis
        t_end = 100
        ts = np.arange(t_begin,t_end,dt_model) 

        # shelf forcing
        Sbottom = 35 # salinity at bottom
        Stop    = 30 # salinity at top
        Tbottom = 4  # temperature at bottom
        Ttop    = 0  # temperature at top
        tw      = 10 # period of oscillation (days)
        zi      = 50+(30/2)*np.sin(2*np.pi*ts/tw) # 'pycnocline' oscillation
        for k in range(len(ts)):
            t_shelf[:,k] = Tbottom-(Tbottom-Ttop)*np.exp(zs/zi[k]); 
            s_shelf[:,k] = Sbottom-(Sbottom-Stop)*np.exp(zs/zi[k]); 

        # in this example there is no subglacial discharge
        tsg = ts
        qsg = 0*tsg
    case 3:
        # time axis
        t_end = 100
        ts = np.arange(t_begin,t_end,dt_model) # time vector for shelf forcing
        
        # shelf forcing
        t_shelf = 3*np.ones([len(zs),len(ts)])
        s_shelf = 34*np.ones([len(zs),len(ts)])
        tsg     = ts    # time vector for subglacial discharge
        qsg     = 0*tsg # subglacial discharge on tsg
        # the icebergs are dealt with in the config file
    case 4:
        # time axis
        t_end   = 2*365
        ts      = np.arange(t_begin,t_end,dt_model) # time vector for shelf forcing

        # shelf forcing
        Sbottom = 35 # salinity at bottom
        Stop    = 30 # salinity at top
        Tbottom = 4 # temperature at bottom
        Ttop    = 0 #temperature at top
        tw      = 10 # period of oscillation (days)
        zi      = 50+(30/2)*np.sin(2*np.pi*ts/tw) # 'pycnocline' oscillation
        for k in range(len(ts)):
            s_shelf[:,k] = Sbottom-(Sbottom-Stop)*np.exp(zs/zi[k]); # shelf salinity
            t_shelf[:,k] = Tbottom-(Tbottom-Ttop)*np.exp(zs/zi[k]); # shelf temperature
            
        # glacier forcing
        tsg = np.arange(t_begin,t_end,dt_model) # create time axis
        qsg = 300*np.exp(-((tsg % 365-200)/30)**2); # subglacial discharge on tsg



n_plumes = 1 # choose number of plumes
i_plumes = [i for i in range(n_plumes)] # plume coordinate

# set up GL depths
gl_depth = 800.0
gl_depth_all_plumes = np.ones(n_plumes)*gl_depth                   # replicate the same GL for all plumes
qsg_all_plumes      = np.repeat(qsg[np.newaxis,:],n_plumes,axis=0) # replicate Qsg for all plumes

#TODO: convert time axes (ts & tsg) to the proper units: starting at 2010-01-01

# set up output datasets
ds_shelf = xr.Dataset(coords=dict(depth=(['depth'],zs),time=(['time'],ts)), \
                      attrs=dict(run_name=f"example{example_number}_subglacial_discharge"))
    
ds_glacier = xr.Dataset(coords=dict(plume=(['plume'],i_plumes),time=(['time'],tsg)), \
                      attrs=dict(run_name=f"example{example_number}_subglacial_discharge"))
    
    
ds_Ts  = xr.DataArray(t_shelf,coords={'depth':zs,'time':ts})
ds_Ss  = xr.DataArray(s_shelf,coords={'depth':zs,'time':ts})
ds_Qsg = xr.DataArray(qsg_all_plumes,coords={'plume':i_plumes,'time':tsg})
ds_Hgl = xr.DataArray(gl_depth_all_plumes,coords={'plume':i_plumes})

ds_shelf['thetao'] = ds_Ts
ds_shelf['so'] = ds_Ss
ds_glacier['Qsg'] = ds_Qsg
ds_glacier['Hgl'] = ds_Hgl

# saving...
ds_shelf.to_netcdf(f"example{example_number}_shelf_forcings.nc")
ds_glacier.to_netcdf(f"example{example_number}_glacier_forcings_np{n_plumes}+'.nc")
