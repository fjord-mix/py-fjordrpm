#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a desired ocean dataset, takes a region subset, and averages it
Saves that as a T and S time series in a NetCDF to be read by FjordRPM

Current set to use Copernicus Marine Services' mean reanalysis product
dataset doi: https://doi.org/10.48670/moi-00024

Created on Tue Dec  3 14:58:12 2024
@author: mmeb1
"""
import xarray     as xr
import datetime   as dt
import juliandate as jd
# %% all needed paths
path_nc = '/Volumes/ice-ocean-group/greenland_common/ocean/'
path_out = '/Users/mmeb1/Library/CloudStorage/OneDrive/data_common/greenland/FjordMIX/processed_data/coupling_tests/'
file_in = 'MP_0p25deg_ocn_grl_julian.nc'
file_out = 'MP_shelf_kang.nc'

# %% coordinates of the box for sampling
lat1 = 67
lat2 = 68
lon1 = -32.4
lon2 = -29.5

# %% opening and subsetting
ds_in = xr.open_dataset(path_nc+file_in)
ds_sel = ds_in.sel(latitude=slice(lat1,lat2),longitude=slice(lon1,lon2)).mean(dim=['latitude','longitude'])

time_gregorian = [dt.datetime(*jd.to_gregorian(t)) for t in ds_sel.time]

ds_out = ds_sel.transpose()
ds_out['time'] = time_gregorian
# %% saving (time,depth) NetCDF
ds_out.to_netcdf(path_out+file_out)
