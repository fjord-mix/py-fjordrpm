#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main wrapper class for running the model. 
Will initialise the model according to the specs in the config file, 
and will run the simulation and save its outputs upon calling "run_fjord()"

Created on Thu Dec  5 10:50:43 2024
@author: mmeb1
"""
import tomli
import fjordrpm_setup   as fst
import fjordrpm_io      as fio
import fjordrpm_physics as fph
from   time           import process_time

class FjordRPM():
    """ 
    Wrapper class for Python's implementation of FjordRPM.
    Loading of config file largely based on the equivalent implementation for LADDIE, credits to Dr. E. (Erwin) Lambert (KNMI)
    """
    def __init__(self,configfile):
        self.startwalltime = process_time()
        with open(configfile,'rb') as f:
            self.config = tomli.load(f)
        fst.init_model(self,configfile)
        
        
    def run_fjord(self):
        self.print2log(f"Starting model run {self.run_name}...")
        self.s = fph.run_model(self.p.to_dict(),self.t,self.f.to_dict(),self.a.to_dict())
        self.print2log("Simulation complete.")
        fio.save_out_nc(self)
        self.print2log("Output saving complete.")
        
    def print2log(self,text):
        """Function to print a line to the log file
        printing of log file largely based on the equivalent implementation for LADDIE, credits to Dr. E. (Erwin) Lambert (KNMI)
        """

        #Ignore if logfile is not defined yet or if we are not meant to write to log
        if hasattr(self,'logfile'):

            #Compute time since start of run
            hours, rem = divmod(process_time()-self.startwalltime, 3600)
            minutes, seconds = divmod(rem, 60)

            with open(self.logfile,"a") as file:
                file.write(f"[{hours:02.0f}:{minutes:02.0f}:{seconds:02.0f}] {text}\n")
            return    

    def throw_error(self,text):
        ''' Essentially a placeholder function for throwing errors
        '''
        self.print2log('ERROR: '+text)
        return
