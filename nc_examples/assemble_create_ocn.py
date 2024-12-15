#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 21:39:07 2024

@author: mmeb1
"""
import xarray as xr
import os


# Assembling all files into one
ds=xr.open_mfdataset("/Volumes/ice-ocean-group/greenland_common/ocean/MP_0.25deg/raw_files/*.nc")
#time_julian=pd.to_datetime(ds.time).to_julian_date() # we need to convert to Julian dates so Matlab can read it
#ds_final_renamed['time']=time_julian
ds_final_renamed=ds.rename_vars({'thetao_mean':'thetao','so_mean':'so'}) # renaming to standardise
ds_final_renamed.to_netcdf("./MP_0p25deg_ocn_grl.nc")

# probably not needed! Remaps to a more manageable grid type for subsetting later
print("Remapping data...")
os.system("cdo remapbil,r1440x720 ./MP_0p25deg_ocn_grl.nc ./MP_0p25deg_ocn_grl_remap.nc")
os.system("rm ./MP_0p25deg_ocn_grl.nc")
os.system("mv /MP_0p25deg_ocn_grl_julian_remap.nc ./MP_0p25deg_ocn_grl.nc")
