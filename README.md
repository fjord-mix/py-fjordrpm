# Python implementation of FjordRPM

The original code written in MATLAB can be found [here](https://github.com/fjord-mix/fjordrpm/tree/main)

# Inputs

How netCDF files should be created:
- **shelf:** variables `thetao` and `so` of dimensions `(depth,time)` with `depth` being positive, and `time` being date in `????`
- **glacier:** `Qsg` with dimensions `(plume,time)` and `Hgl` with dimenions `(plume)`, `plume` being which plume (starting from zero) and `time` being date in `????`

All other parameters can be set up in the config file. The folder `examples` contains config files equivalent to the 4 examples in the MATLAB repository. The folder `nc_examples` contains example scripts for creating the forcings (both idealised and from ocean reanalysis/runoff datasets).

# Outputs
Saves a NetCDF file with all variables. Plume-related quantities have dimensions `(plume,depth,time)`, whereas all others have dimensions `(depth,time)`. Certain time-independent or depth-indepdendent variables will not have the respective dimension.
