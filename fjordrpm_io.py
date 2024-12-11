#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set of functions for handling input/output of the model, including parsing the config file

Several functions were adapted from LADDIE's source code, 
copyright by Dr Erwin Lambert (KNMI) distributed under a BSD-3 Clause license

Created on Thu Dec  5 10:52:18 2024
@author: mmeb1
"""
import os,sys
import datetime as dt
import pandas   as pd
import fjordrpm_utils as fut

def tryread(object,category,parameter,reqtype,valid=None,allowconversion=True,checkfile=False,checkdir=False,default=None):
    """
    Function to read values from config-file, check type and values, and aborting or defaulting if missing
    This function was copied from LADDIE's source code, 
    copyright by Dr Erwin Lambert (KNMI) distributed under a BSD-3 Clause license
    """
    

    #Make sure integers 0 or 1 are not interpreted as boolean
    if reqtype==bool:
        allowconversion=False

    #Check whether input parameter exists
    try:
        out = object.config[category][parameter]
    except:
        if default == None:
            print(f"INPUT ERROR: missing input parameter '{parameter}' in [{category}]. Please add to config-file of type {reqtype.__name__}")
            sys.exit()
        else:
            if object.newdir:
                object.print2log(f"Note: missing input parameter '{parameter}' in [{category}], using default value {default}")
            out = default

    #Check whether input parameter is of the correct type. Second criterium because bool is treated as subtype of integer in python
    if (isinstance(out,reqtype) == False or (isinstance(out,bool) and reqtype == int)):
        if allowconversion:
            try:
                #Convert to required type, for example float to int or vice versa
                out2 = reqtype(out)
                if object.newdir:
                    object.print2log(f"Note: changing input parameter '{parameter}' from {type(out).__name__} to {reqtype.__name__}")
                out = out2
            except:
                if default == None:
                    print(f"INPUT ERROR: input parameter '{parameter}' in [{category}] is of wrong type. Is {type(out).__name__}, should be {reqtype.__name__}")
                    sys.exit()
                else:
                    print(f"WARNING: wrong type '{parameter}' in [{category}], using default value {default}")
                    out = default
        else:
            if default == None:
                print(f"INPUT ERROR: input parameter '{parameter}' in [{category}] is of wrong type. Is {type(out).__name__}, should be {reqtype.__name__}")
                sys.exit()            
            else:
                print(f"WARNING: wrong type '{parameter}' in [{category}]. Is {type(out).__name__}, should be {reqtype.__name__}, using default value {default}")
                out = default          

    #Check whether value of input is valid
    if valid != None:
        if isinstance(valid,list):
            if out not in valid:
                if default == None:
                    print(f"INPUT ERROR: invalid value for '{parameter}' in [{category}]; choose from {valid}")
                    sys.exit()
                else:
                    print(f"WARNING: invalid value '{parameter}' in [{category}], using default value {default}")
                    out = default
        if isinstance(valid,tuple):
            if out < valid[0]:
                if default == None:
                    print(f"INPUT ERROR: invalid value for '{parameter}' in [{category}]; should be >= {valid[0]} ")
                    sys.exit()
                else:
                    print(f"WARNING: invalid value '{parameter}' in [{category}]; should be >= {valid[0]}, using default value {default}")
                    out = default
            if out > valid[1]:
                if default == None:
                    print(f"INPUT ERROR: invalid value for '{parameter}' in [{category}]; should be <= {valid[1]} ")
                    sys.exit()
                else:
                    print(f"WARNING: invalid value '{parameter}' in [{category}]; should be <= {valid[1]}, using default value {default}")
                    out = default

    #Check whether file exists
    if checkfile:
        if os.path.isfile(out) == False:
            print(f"INPUT ERROR: non-existing file for '{parameter}' in [{category}]; check filename")
            sys.exit()
        if out[-3:] != ".nc":
            print(f"INPUT ERROR: file '{parameter}' in [{category}] must be '.nc'; check filename")
            sys.exit()

    #Check whether directory exists
    if checkdir:
        if os.path.isdir(out) == False:
            try:
                os.mkdir(out)
                print('WARNING: making a new results directory')
            except:
                print(f"INPUT ERROR: could not create directory '{parameter}' in [{category}]; check directory name")
                sys.exit()

    return out


def create_out_dir(object,configfile):
    """Create run directory and logfile
    This function was copied from LADDIE's source code, 
    copyright by Dr Erwin Lambert (KNMI) distributed under a BSD-3 Clause license
    """

    #First assume new dir is created
    object.newdir = True

    #Read run name
    object.name = tryread(object,"Run","run_name",str)
    #Read directory to store output
    object.resultdir = tryread(object,"Outputs","output_dir",str,checkdir=True)
    #Read flag to force new folder
    object.forcenewdir = tryread(object,"Outputs","force_new_dir",bool,default=True)
    #Read logfile
    object.logfilename = tryread(object,"Outputs","log_file",str,default="log.txt")

    #Derive desired run directory:
    object.rundir = os.path.join(object.resultdir,object.name)
    if object.forcenewdir or os.path.isdir(object.rundir) == False:
        try:
            #Create rundirectory
            os.mkdir(object.rundir)
        except:
            try:
                #Create rundirectory with current date
                object.rundir = os.path.join(object.resultdir,dt.datetime.today().strftime(f"{object.name}_%Y-%m-%d"))
                os.mkdir(object.rundir)
            except:
                for n in range(100):
                    try:
                        #Create rundirectory with current date and incremental number. Give up after 100 tries
                        object.rundir = os.path.join(object.resultdir,dt.datetime.today().strftime(f"{object.name}_%Y-%m-%d_{n}"))
                        os.mkdir(object.rundir)
                        break
                    except:
                        continue
    else:
        #No new directory is created, but using existing directory to continue run
        object.newdir = False

    #Create log file
    object.logfile = os.path.join(object.rundir,object.logfilename)

    #Copy config file to run directory
    os.system(f"cp {configfile} {object.rundir}")

    if object.newdir:
        object.print2log(f"Created new run directory {object.rundir}")
    else:
        object.print2log(f"Continuing run in existing directory {object.rundir}")

    return

def parse_config(object):
    """Inherit all config keys and check whether they are of the correct form / type
    This function was adapted from LADDIE's source code, 
    copyright by Dr Erwin Lambert (KNMI) distributed under a BSD-3 Clause license
    
    """
    #object.save_log    = tryread(object,"Outputs","save_log",bool,default=True)

    object.print2log("========= Starting to read config file =============")

    #Run
    object.run_name = tryread(object,"Run","run_name",str)

    #Time
    object.date_begin = tryread(object,"Time","time_begin",dt.date)
    object.time_begin = pd.to_datetime(object.date_begin).to_julian_date()
    object.date_end   = tryread(object,"Time","time_end",dt.date)
    object.time_end   = pd.to_datetime(object.date_end).to_julian_date()
    object.dt_model   = tryread(object,"Time","dt_model",float)
    object.dt_plume   = tryread(object,"Time","dt_plume",float)


    #Forcing
    object.shelf_forcing_file   = tryread(object,"Forcings","shelf_forcing_file",str,checkfile=False,default='')
    object.glacier_forcing_file = tryread(object,"Forcings","glacier_forcing_file",str,checkfile=False,default='')

    object.iceberg_profile = tryread(object,"Forcings","iceberg_profile",str,["exponential1","exponential2"])
    match object.iceberg_profile:
        case 'exponential1':
            #TODO: implement symbolic maths for this?  
            object.iceberg_fun = fut.iceberg_fun_exponential1

        case 'exponential2':
            #TODO: implement symbolic maths for this?
            object.iceberg_fun = fut.iceberg_fun_exponential2
        case _:
            #TODO: implement symbolic maths for this?
            object.iceberg_fun = fut.iceberg_fun_exponential1

    #Initialisation    
    object.start_mode   = tryread(object,"Initialisation","start_mode",str,['init_from_t0','init_from_mean','restart_from_run'])
    object.restart_file = tryread(object,"Initialisation","restart_file",str,checkfile=False,default='')

    #Geometry
    object.L          = tryread(object,"Geometry","fjord_length",float,(0,1e20))
    object.W          = tryread(object,"Geometry","fjord_width",float,(0,1e20))
    object.H          = tryread(object,"Geometry","fjord_depth",float,(0,1e20))
    object.Hgl        = tryread(object,"Geometry","fjord_gl_depth",float,(0,1e20))
    object.Hsill      = tryread(object,"Geometry","fjord_sill_depth",float,(0,1e20))
    object.hypsometry = tryread(object,"Geometry","type_hypsometry",str,['cuboid','formula','read_from_file'])
    object.N          = tryread(object,"Geometry","fjord_n_layers",int,(0,1e20))
    match object.hypsometry:
        case 'cuboid':
            pass # nothing to do
        case 'idealised':
            object.hypsometry_fun = fut.hypsometry_idealised
            #TODO: implement idealised function?
            pass
        case 'read_from_file':
            object.hypsometry_file = tryread(object,"Geometry","hypsometry_file",str,checkfile=False,default='')
            pass
        case _:
            pass # same as cuboid
    
    #Parameters
    object.A0     = tryread(object,"Parameters","iceberg_area_A0",float,(0,1e20),default=3e8)
    object.P0     = tryread(object,"Parameters","plume_width_P0",float,(0,1e20),default=250)
    object.C0     = tryread(object,"Parameters","shelf_exchange_C0",float,(0,1e20),default=1e4)
    object.K0     = tryread(object,"Parameters","vertical_mixing_K0",float,(0,1e20),default=5e-3)
    object.Kb     = tryread(object,"Parameters","background_mixing_Kb",float,(0,1e20),default=1e-5)
    object.Ri0    = tryread(object,"Parameters","richardson_mixing_Ri0",float,(0,1e20),default=700)
    object.M0     = tryread(object,"Parameters","iceberg_melt_M0",float,(0,1e20),default=5e-7)
    object.U0     = tryread(object,"Parameters","iceberg_upwelling_scale_U0",float,(0,1),default=1)
    object.nu     = tryread(object,"Parameters","iceberg_decay_scale_nu0",float,(0,1e20),default=25)   
    

    # Physical constants
    object.g      = tryread(object,"Constants","g",float,(0,1e20),default=9.81)
    object.betaS  = tryread(object,"Constants","betaS",float,(0,1e20),default=7.86e-4)
    object.betaT  = tryread(object,"Constants","betaT",float,(0,1e20),default=3.87e-5)
    object.l      = tryread(object,"Constants","l",float,(0,1e20),default=3.34e5)
    object.cp     = tryread(object,"Constants","cw",float,(0,1e20),default=3974)
    object.ci     = tryread(object,"Constants","ci",float,(0,1e20),default=2009)
    object.l1     = tryread(object,"Constants","l1",float,(-1e20,0),default=-5.73e-2)
    object.l2     = tryread(object,"Constants","l2",float,(0,1e20),default=8.32e-2)
    object.l3     = tryread(object,"Constants","l3",float,(-1e20,0),default=-7.61e-4) # is this meant to be positive or negative?
    object.GT     = tryread(object,"Constants","GT",float,(0,1e20),default=2.2e-2)
    object.GS     = tryread(object,"Constants","GS",float,(0,1e20),default=6.2e-4)
    object.Cd     = tryread(object,"Constants","Cd",float,(0,1e20),default=2.5e-3)
    object.Ti     = tryread(object,"Constants","Ti",float,(-1e20,0),default=-10)
    object.alphai = tryread(object,"Constants","alphai",float,(0,1),default=0.1)
    object.alphap = tryread(object,"Constants","alphap",float,(0,1),default=0.1)
    object.sid    = tryread(object,"Constants","sid",int,(0,1e20),default=86400)
    
    #Outputs
    #object.output_dir = tryread(object,"Outputs","output_dir",str)


    object.print2log("============= Finished reading config. All input correct =======")
    object.print2log("================================================================")
    return

def save_out_nc(object): #TODO: implement saving routine (after model physics is translated?) 
    object.print2log(f"Saving outputs in {object.resultdir}...")  
    #create_DataFrame
    #add_outputs
    #add_forcings
    #add_geometry
    #convert_time_axis_from_julian
    #write_ds_to_disk
    return