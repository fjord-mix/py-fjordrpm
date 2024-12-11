#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for triggering a single model run
This function was largely adapted from LADDIE's source code, 
copyright by Dr Erwin Lambert (KNMI) distributed under a BSD-3 Clause license

Created on Fri Dec  6 12:51:30 2024
@author: mmeb1
"""
import sys
# =============================================================================
# sys.path.append('src/')
# =============================================================================
from fjordrpm import FjordRPM

#Check whether config file is specified
#Nargs = len(sys.argv)
#assert Nargs > 1, f"Need to specify config file"
#assert Nargs < 3, f"Too many arguments"
#assert sys.argv[1][-4:] == 'toml' , f"config file should be .toml"

#fjord_run = Feshie(sys.argv[1])
fjord_run = FjordRPM('config_test_kang.toml')
fjord_run.run_fjord()