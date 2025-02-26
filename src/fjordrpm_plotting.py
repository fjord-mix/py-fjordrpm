#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 18:44:07 2024

@author: mmeb1
"""
import xarray as xr

path_in = '/Users/mmeb1/fjordrpm_coupling/outputs/example_restart_from3/'
file_in = 'outputs_example_restart_from3.nc'


ds = xr.open_dataset(path_in+file_in)

# OK, depth needs to be negative...
