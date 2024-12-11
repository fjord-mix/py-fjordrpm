#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of functions to setup the model run: 
    1. reads and parses the config file
    2. creates the output directory
    3. creates the main structures/objects: (t)ime,(p)parameters,(f)forcings, initi(a)l conditions

Created on Wed Dec  4 16:16:58 2024
@author: mmeb1
"""
import numpy     as np
import fjordrpm_io as fio
from fjordrpm_classes import Parameters, Forcings, InitialState


def create_time_axis(object):
    t = np.arange(object.time_begin,object.time_end,object.dt_model)
    return t

def load_inputs(object):    
    object.t = create_time_axis(object)
    object.print2log("Julian time axis created.")

    object.p = Parameters(object)
    object.print2log("Model parameters successfully set up.")
    
    object.f = Forcings(object)
    if hasattr(object.f,'Hgl'):
        object.p.set_gl(object.f.Hgl) # if .nc file has Hgl, it will supersede the one in the config file
    object.print2log("Model forcings successfully set up.")
    object.a = InitialState(object)
    if hasattr(object,'hypsometry_fun') or hasattr(object,'hypsometry_file'):
        object.p.set_hypsometry(object,object.a.z)

        
    object.print2log("Initial conditions successfully set up.")
    return

def init_model(object,configfile):
    object.print2log("============= Initialising model =======")
    fio.parse_config(object)
    fio.create_out_dir(object,configfile)
    load_inputs(object)
    return

    
    
    
