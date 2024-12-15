#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes for Parameters (p), Forcings (f), and InitialState (a) of the model
They will ensure all forcings and options in the config file will be converted to the model's input style

Created on Thu Dec  5 10:56:11 2024
@author: mmeb1
"""
import xarray as xr
import numpy  as np
import pandas as pd
import fjordrpm_utils as fut
from scipy.interpolate import interp1d

class Parameters:
    
    def __init__(self,config):
        self.A0     = config.A0
        self.wp     = np.array([config.wp])
        self.C0     = config.C0
        self.K0     = config.K0
        self.Kb     = config.Kb
        self.Ri0    = config.Ri0
        self.M0     = config.M0
        self.U0     = config.U0
        self.nu     = config.nu
        self.g      = config.g
        self.betaS  = config.betaS
        self.betaT  = config.betaT
        self.l      = config.l
        self.cw     = config.cw
        self.ci     = config.ci
        self.l1     = config.l1
        self.l2     = config.l2
        self.l3     = config.l3
        self.GT     = config.GT
        self.GS     = config.GS
        self.Cd     = config.Cd
        self.Ti     = config.Ti
        self.alphai = config.alphai
        self.alphap = config.alphap
        self.sid    = config.sid
        
        self.t_save          = np.arange(config.time_begin,config.time_end,config.dt_out)
        self.run_plume_every = np.floor(config.dt_model/config.dt_plume)
        
        # geometry attributes
        self.L      = config.L
        self.H      = config.H
        self.Hgl    = config.Hgl
        self.Hsill  = config.Hsill
        self.W      = config.W    
        self.N      = config.N        
        self.sill   = config.sill
        
        if self.Hsill >= self.H:
            self.sill = 0
      
        return
    
    # Hsill might need updates if glacier NetCDF has n_plumes grounding-line depths
    def set_gl(self,gl_depths):
        self.Hgl = gl_depths
        if len(gl_depths) > 1 and len(self.wp) ==1:
            wp = self.wp
            self.wp = np.array([wp for _ in range(len(gl_depths))])

    # in case we need to update the plume widths from say a NetCDF
    def set_wp_array(self,wp_array):
        self.wp = wp_array
        return

    def set_hypsometry(self,object,z):
        if hasattr(object,'hypsometry_fun'):
            self.W = object.config.hypsometry_fun(z)
        elif hasattr(object,'hypsometry_file'):
            ds_hyp = xr.open_dataset(object.object.hypsometry_file)
            self.W = np.interp(z,ds_hyp.depth,ds_hyp.width)
        return
    
    
    def to_dict(self):
        '''returns all relevant attributes as a single dictionary'''
        p = {}
        p['A0']  = self.A0
        p['C0']  = self.C0
        p['wp']  = self.wp
        p['K0']  = self.K0
        p['Kb']  = self.Kb
        p['Ri0'] = self.Ri0
        p['M0']  = self.M0
        
        p['U0']     = self.U0    
        p['nu']     = self.nu    
        p['g']      = self.g     
        p['betaS']  = self.betaS 
        p['betaT']  = self.betaT 
        p['l']      = self.l     
        p['cw']     = self.cw
        p['ci']     = self.ci    
        p['l1']     = self.l1    
        p['l2']     = self.l2    
        p['l3']     = self.l3    
        p['GT']     = self.GT    
        p['GS']     = self.GS    
        p['Cd']     = self.Cd    
        p['Ti']     = self.Ti    
        p['alphai'] = self.alphai
        p['alphap'] = self.alphap
        p['sid']    = self.sid           
        p['W']      = self.W     
        p['L']      = self.L     
        p['H']      = self.H     
        p['Hgl']    = self.Hgl   
        p['Hsill']  = self.Hsill 
        p['N']      = self.N     
        p['sill']   = self.sill
        
        p['t_save']          = self.t_save
        p['run_plume_every'] = self.run_plume_every
      
        return p

class Forcings:
    
    def get_shelf_forcing(self,config):
        ds = xr.open_dataset(config.shelf_forcing_file)
        # make changes to formatting to ensure the proper Ts and Ss structure
        #isinstance(time_begin,pd.Timestamp)
        self.ts = pd.to_datetime(ds.time).to_julian_date()
        self.Ts = ds.thetao
        self.Ts['time'] = self.ts
        self.Ss = ds.so
        self.Ss['time'] = self.ts
        self.zs = ds.depth #TODO: should depths be positive or negative? bottom-up or top-down? or does it matter at all?
        
        # checking appropriate dimensions
        if (self.Ts.shape[0] != self.zs.shape[0]) or (self.Ts.shape[1] != self.ts.shape[0]):
            if (self.Ts.shape[0] == self.ts.shape[0]):
                self.Ts = self.Ts.T
            else:
                #TODO: throw error
                return
            if (self.Ts.shape[1] == self.zs.shape[0]):
                self.Ts = self.Ts.T
            else:
                #TODO: throw error
                return
        if (self.Ss.shape[0] != self.zs.shape[0]) or (self.Ss.shape[1] != self.ts.shape[0]):
            if (self.Ss.shape[0] == self.ts.shape[0]):
                self.Ss = self.Ss.T
            else:
                #TODO: throw error
                return
            if (self.Ss.shape[1] == self.zs.shape[0]):
                self.Ss = self.Ss.T
            else:
                #TODO: throw error
                return    
            
            
        return 
    
    def get_glacier_forcing(self,config):
        ds = xr.open_dataset(config.glacier_forcing_file)
        # make changes to formatting to ensure the proper Qsg and Hgl structure
        self.Qsg = ds.Qsg
        self.tsg = pd.to_datetime(ds.time).to_julian_date()
        try:
            self.Hgl = np.abs(ds.Hgl) # ensures GL depths are positive
        except:
            #TODO: throw warning that will be using Hgl from config_file
            pass
        #check dimensions
        if (self.Qsg.shape[0] != self.Hgl.shape[0]) or (self.Qsg.shape[1] != self.tsg.shape[0]):
            if (self.Qsg.shape[1] == self.Hgl.shape[0]):
                self.Qsg = self.Qsg.T
            else:
                #TODO: throw error
                return
            if (self.Qsg.shape[0] == self.tsg.shape[0]):
                self.Qsg = self.Qsg.T
            else:
                #TODO: throw error
                return
        return
    

    def __init__(self,config):
        self.get_shelf_forcing(config)
        self.get_glacier_forcing(config)
        return

    def to_dict(self):
        '''returns all relevant attributes as a single dictionary'''
        f = {}
        f['Ts'] = self.Ts
        f['Ss'] = self.Ss
        f['zs'] = self.zs
        f['ts'] = self.ts
        f['Qsg'] = self.Qsg
        f['tsg'] = self.tsg
        return f


class InitialState:
    def __init__(self,object):
        match object.start_mode:
            case 'init_from_t0':
                self.init_from_t0(object)
                object.print2log(f"Initial T,S conditions set as closest shelf conditions to {object.date_begin}")
            case 'init_from_mean':
                self.init_from_avg(object)
                object.print2log(f"Initial T,S conditions set as mean shelf conditions between {object.date_begin} and {object.date_end}")
            case 'restart_from_run':
                self.init_from_restart(object.config)
                object.print2log("Initial T,S conditions set from a previous simulation:")
                object.print2log(f"{object.restart_file}")
            case _:
                self.init_from_t0(object)    
                object.print2log(f"Initial T,S conditions set as closest shelf conditions to {object.date_begin}")

        # Getting depths intervals for the fjord
        self.H0 = np.ones([object.p.N,]) * object.p.H/object.p.N
        intervals = np.cumsum(self.H0)
        intervals = np.insert(intervals,0,0)
        intervals = (intervals[1:]+intervals[0:-1])/2
        self.z = intervals

        self.T0 = self.T0.ffill(dim='depth')
        self.S0 = self.S0.ffill(dim='depth')
        # interpolate initial conditions to the fjord depths
        self.T0 = np.interp(self.z, object.f.zs, self.T0)
        self.S0 = np.interp(self.z, object.f.zs, self.S0)

        self.get_iceberg_forcing(object)
        object.print2log("Iceberg area profile set.")
        return
    
    def init_from_t0(self,object):
        self.T0 = object.f.Ts.sel(time=fut.nearest(object.f.ts,object.time_begin))
        self.S0 = object.f.Ss.sel(time=fut.nearest(object.f.ts,object.time_begin))
        return
    
    def init_from_avg(self,object):
        self.T0 = object.f.Ts.sel(time=slice(fut.nearest(object.f.ts,object.time_begin),fut.nearest(object.f.ts,object.time_end))).mean(dim='time',skipna=True)
        self.S0 = object.f.Ss.sel(time=slice(fut.nearest(object.f.ts,object.time_begin),fut.nearest(object.f.ts,object.time_end))).mean(dim='time',skipna=True)
        return
    
    def init_from_restart(self,config): #TODO:implement restart option in config_file
        pass
    
    def get_iceberg_forcing(self,object):
        if object.config.iceberg_profile == 'exponential1':
            iceprofile = object.iceberg_fun(object.nu,object.Hgl,-self.z)
            self.I0 = (object.A0/np.sum(iceprofile))*iceprofile
        else:
            self.I0 = object.A0 * object.iceberg_fun(object.nu,object.Hgl,-self.z)
        return

    def to_dict(self):
        '''returns all relevant attributes as a single dictionary'''
        a = {}
        a['T0'] = self.T0
        a['S0'] = self.S0
        a['I0'] = self.I0
        a['H0'] = self.H0
        return a
        
        
        
        
