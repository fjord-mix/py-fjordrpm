#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 18:44:07 2024

@author: mmeb1
"""
import xarray as xr

path_in = '/Users/mmeb1/fjordrpm_coupling/test_inputs/example3_icebergs/'
file_in = 'outputs_example3_icebergs.nc'


ds = xr.open_dataset(path_in+file_in)

# OK, depth needs to be negative...