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
sys.path.append('src/')
from fjordrpm import FjordRPM

#Check whether config file is specified
Nargs = len(sys.argv)
assert Nargs > 1, "Need to specify config file"
assert Nargs < 3, "Too many arguments"
assert sys.argv[1][-4:] == 'toml' , "config file should be .toml"

fjord_run = FjordRPM(sys.argv[1])
#fjord_run = FjordRPM('examples/config_example1.toml')
fjord_run.run_fjord()

# Use the block below if profiling the code
# =============================================================================
#import cProfile
#import pstats
# cProfile.run('fjord_run.run_fjord()','./profiling_'+sys.argv[1][:-4])
# p = pstats.Stats('./profiling_'+sys.argv[1][:-4])
# p.strip_dirs()
# p.sort_stats(pstats.SortKey.CUMULATIVE)
# p.print_stats(10)
# =============================================================================
