# Python implementation of FjordRPM

This is the Python implementation of the fjord reduced physics model, FjordRPM. 
The original code written in MATLAB can be found [here](https://github.com/fjord-mix/fjordrpm/tree/main). Please refer to that repository and the model description paper for most of the information about FjordRPM itself.
**Disclaimer:** This is more of a short side project than an official release! Although a lot of pre- and post-processing code has been written from scratch (and some taken from [LADDIE](https://github.com/erwinlambert/laddie), which is FjordRPM's distant Antarctic cousin), the model core and plume model were initially translated using ChatGPT to save time, but the code has been extensively reviewed afterwards. Luckily for us, ChatGPT is far from being perfect and we coders are still required. Since this is not the official release, the code has only been thoroughly tested and compared against the examples in the MATLAB repository. Do not hesitate, however, to get in touch should you spot any bugs or have any suggestions for improvement!

The main differences between pyFjordRPM and FjordRPM are summarised below, but mainly revolve around how inputs and outputs work.

# Usage:

Once FjordRPM is installed (see [Installation](#Installation)), simply type `python run_fjordrpm.py <path_to_config_file>` in the terminal

# Inputs

FjordRPM will run most model parameters from the specified config file, including paths to the shelf forcing and glacier forcing NetCDF files.

How input netCDF files should be created:
- **shelf:** variables `thetao` and `so` for shelf temperature and salinity, respectively, of dimensions `(depth,time)` with `depth` being positive, and `time` being an array of `datetime` objects
- **glacier:** `Qsg` with dimensions `(plume,time)` and `Hgl` with dimenions `(plume)`, `plume` being which plume (starting from zero) and `time` being an array of `datetime` objects

All other parameters can be set up in the config file. The folder `examples` contains config files equivalent to the 4 examples in the MATLAB repository. The folder `nc_examples` contains example scripts for creating the forcings (both idealised and from ocean reanalysis/runoff datasets).

# Outputs
Saves a NetCDF file with all variables. Plume-related quantities have dimensions `(plume,depth,time)`, whereas all others have dimensions `(depth,time)`. Certain time-independent or depth-indepdendent variables will not have the respective dimension.
Outputted files can be used to take the end of a previous run as initial conditions by setting up the `[Initialisation]` part of the config file as:
```
start_mode   = 'restart_from_run'
restart_file = 'path/to/output_file_from_previous_run.nc' 
```

# Installation

1. Download or clone the repository by using `git clone https://github.com/fjord-mix/py-fjordrpm`
2. Set up a virtual environment by using `conda env create -f fjordrpm.yaml`. This environment should contain all needed python packages.