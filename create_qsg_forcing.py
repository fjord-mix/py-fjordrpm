#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in .csv files from Karlsson et al. (2023): Greenland Freshwater Flux on Glacier-basin Scales
dataset doi: https://doi.org/10.22008/FK2/BOVBVR
Saves that as a .nc file to be read by FjordRPM

Created on Tue Dec  3 15:12:27 2024
@author: mmeb1
"""
import pandas as pd
import xarray as xr
import glob

# %% Paths
path_files = '/Volumes/ice-ocean-group/greenland_common/runoff/Karlsson2023/'
path_out = '/Users/mmeb1/Library/CloudStorage/OneDrive/data_common/greenland/FjordMIX/processed_data/coupling_tests/'
file_out = 'Kangerlussuaq_Karlsson2023.nc'

# %% the actual relevant stuff we would need as input
gate_IDs = [188,190,192,194]
Hgl      = [100, 200, 300, 400]

# %% getting the file names
gl_depths = xr.DataArray(Hgl,dims='plume',coords={'plume':range(len(gate_IDs))})
files_in = [glob.glob(path_files+'*_D'+str(gate_id)+'.csv') for gate_id in gate_IDs]

# %% creating NetCDF from all flux gate CSVs
is_first = True
for file_in in files_in:
    ds_in = pd.read_csv(file_in[0],skiprows=2)
    ds_in.rename(columns={'Date_YYYY-MM':'time'}, inplace=True)
    ds_in.index = pd.to_datetime(ds_in['time'])
    
    
    days_in_month = pd.to_datetime(ds_in.time).dt.daysinmonth # we will need that to convert monthly to daily values
    days_in_month = days_in_month.asfreq(freq='D',method='ffill') # we use forward fill (ffill) since to_datetime converts YYYY-MM as the first day of the month
    
    
    ds_in['Qsg'] = ds_in['SurfaceMelt']+ds_in['BasalMelt'] # this is what we actually want
    ds_qsg = ds_in[['Qsg']] # get rid of 'datetime' column that we do not need anymore (only as index)
    
    
    ds_daily = ds_qsg.asfreq(freq='D',method='ffill')
    ds_daily = ds_daily.divide(days_in_month,axis=0)
    ds_daily = ds_daily.to_xarray().expand_dims(dim='plume',axis=0)        
    
    if is_first:
        ds_out = ds_daily
        is_first=False
    else:
        ds_out = xr.concat([ds_out,ds_daily],dim='plume')

# %% Adding grounding-line depths for each gate/plume        
ds_out['Hgl'] = gl_depths

# %% saving output        
ds_out.to_netcdf(path_out+file_out)


