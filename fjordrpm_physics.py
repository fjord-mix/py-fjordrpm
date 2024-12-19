#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:05:45 2024

@author: original code from Slater et al. (2025; GMD), translated to Python by ChatGPT and reviewed by Martim Mas e Braga
"""
import numpy as np
from scipy.interpolate import interp1d
import fjordrpm_fluxes as ffl
import fjordrpm_utils  as fut


def check_inputs(p, t, f, a):
    """
    CHECK_INPUTS Check for errors in model inputs.
    status = CHECK_INPUTS(p, t, f, a) checks the model inputs.
    """
    status = 0

    # Check provided depths are positive
    if any(depth < 0 for depth in [p['Hsill'], p['H']]) or any(p['Hgl'] < 0):
        status = 1
        raise ValueError('p.H, p.Hgl and p.Hsill must be positive')

    # Check dimensionality of initial conditions
    if not (a['H0'].shape == a['S0'].shape == a['T0'].shape == a['I0'].shape == (p['N'],)):
        status = 1
        raise ValueError('Initial conditions (a.H0, a.S0, a.T0, and a.I0) must have dimensions p.N x 1')


    # Check sum of layer thicknesses is equal to fjord depth
    if abs(np.sum(a['H0']) - p['H']) > 1e-10:
        status = 1
        raise ValueError('Layer thicknesses (a.H0) must sum to fjord depth (p.H)')

    # Check dimensionality of shelf forcing
    nz = len(f['zs'])
    nt = len(f['ts'])


    if f['Ss'].shape != f['Ts'].shape or f['Ss'].shape != (nz, nt):
        status = 1
        raise ValueError('f.Ss and f.Ts must have dimensions length(f.zs) x length(f.ts)')

    # Check dimensionality of discharge forcing
    nt = len(f['tsg'])


    if f['Qsg'].shape[1] != nt:
        status = 1
        raise ValueError('Second dimension of f.Qsg must have length nt')

    # Check number of plumes
    if f['Qsg'].shape[0] != len(p['wp']) or f['Qsg'].shape[0] != len(p['Hgl']):
        status = 1
        raise ValueError('Check num_plumes==f.Qsg.shape[0]==len(p.wp)==len(p.Hgl)')

    # Uncomment the below section if you want to check if t_save is a subset of t
    if 't_save' not in p:
        p['t_save'] = t
    if not np.all(np.isin(p['t_save'], t)):
        status = 1
        raise ValueError("The values where the solution is saved, p.t_save, must be a subset of the values where the solution is computed, t.")

    return status


def initialise_variables(p, t, f, a):
    """
    INITIALISE_VARIABLES Initialise model variables with initial conditions.
    s = INITIALISE_VARIABLES(p, t, f, a) initialises the model variables
    using parameters structure p, time vector t, and initial conditions
    structure a, and returns solution structure s.
    """
    s = {}

    # Set the timestep from the input time field
    s['dt'] = t[1:] - t[:-1]

    # Initialise solution structure fields
    # Fields with dimensions p.N x 1
    s['H'], s['V'], s['I'] = np.zeros((p['N'], 1)), np.zeros((p['N'], 1)), np.zeros((p['N'], 1))

    # Fields with dimensions p.N x length(t)
    s['T'], s['S'], \
    s['QVs'], s['QTs'], s['QSs'], s['Ts'], s['Ss'], s['phi'], \
    s['QVk'], s['QTk'], s['QSk'], \
    s['QVi'], s['QTi'], s['QSi'], s['QMi'], \
    s['QVv'], s['QTv'], s['QSv'] = [np.zeros((p['N'], len(t))) for _ in range(18)]

    # Fields with dimensions num plumes x p.N x length(t)
    num_plumes = len(p['wp'])
    s['QVp'], s['QTp'], s['QSp'], s['QMp'], s['QEp'] = \
        [np.zeros((num_plumes, p['N'], len(t))) for _ in range(5)]

    # Fields with dimensions num plumes x length(t)
    s['knb'] = np.zeros((num_plumes, len(t)),dtype=int) * np.nan
    s['Qsg'] = np.zeros((num_plumes, len(t)))

    # Initialise layer depths 
    # first deal with case where sill is so shallow or so deep that it would
    # result in less than half a layer at top or bottom, by tweaking the
    # sill depth itself
    if p['Hsill'] < p['H'] / p['N']:  # Avoid very thin layers at top
        p['Hsill'] = p['H'] / p['N']
    elif p['Hsill'] >= p['H'] - 0.5 * p['H'] / p['N']:  # Very deep sill, round to no sill
        p['Hsill'] = p['H']
        p['sill'] = 0
    elif p['Hsill'] >= p['H'] - p['H'] / p['N']:  # Avoid very thin layers at bottom
        p['Hsill'] = p['H'] - p['H'] / p['N']

    # Make sill depth coincide with layer boundary
    if p['sill'] == 1:
        # Redistribute the layers above and below the sill depth
        Nabove = round((p['Hsill'] / p['H']) * p['N'])
        Nbelow = p['N'] - Nabove
        s['H'] = np.concatenate([(p['Hsill'] / Nabove) * np.ones(Nabove),
                                 ((p['H'] - p['Hsill']) / Nbelow) * np.ones(Nbelow)])
        s['ksill'] = Nabove
    else:
        s['H'] = a['H0']
        s['ksill'] = p['N']-1

    # Initialise layer volumes
    s['V'] = s['H'] * p['W'] * p['L']

    # Get forcings on model layers and at model time steps
    s['Ts'], s['Ss'], s['Qsg'] = fut.bin_forcings(f, s['H'], t)

    # Set any discharge values less than 1e-3 to 0
    s['Qsg'][s['Qsg'] < 1e-3] = 0

    # Find the layer with the grounding line and store its index
    ints = np.cumsum(s['H'])
    s['kgl'] = np.zeros(len(p['Hgl']), dtype=int)
    for j in range(len(p['Hgl'])):
        s['kgl'][j] = np.where(ints == fut.nearest(ints,p['Hgl'][j]))[0][0]

    # Redistribute initial conditions according to the new layer boundaries
    ints_old = np.concatenate(([0], np.cumsum(a['H0'])))
    centres_old = 0.5 * (ints_old[:-1] + ints_old[1:])
    ints_new = np.concatenate(([0], np.cumsum(s['H'])))
    centres_new = 0.5 * (ints_new[:-1] + ints_new[1:])

    # Interpolate initial conditions (T, S, I) to the new layer centers
    s['T'][:, 0] = np.interp(centres_new, centres_old, a['T0'])
    s['S'][:, 0] = np.interp(centres_new, centres_old, a['S0'])
    s['I'] = np.interp(centres_new, centres_old, a['I0'])

    return s


def homogenise_unstable_layers(i, p, s):
    """
    HOMOGENISE_UNSTABLE_LAYERS Homogenizes two unstably stratified layers.
    s = HOMOGENISE_UNSTABLE_LAYERS(i, p, s) takes input timestep i,
    parameters structure p, and solution structure s and computes the
    buoyancy jump between layers in the solution for that timestep,
    homogenizes any layers that are unstably stratified and returns an
    updated solution structure s.

    Parameters:
    - i: Current timestep index
    - p: Dictionary containing model parameters (e.g., g, betaS, betaT)
    - s: Dictionary containing model solution at current timestep (e.g., V, T, S)

    Returns:
    - Updated solution dictionary s with homogenized temperature and salinity
    """
    # Required variables
    V = s['V']
    T = s['T'][:, i]  # Temperature for current timestep
    S = s['S'][:, i]  # Salinity for current timestep

    # Compute the buoyancy jump between layers
    B = p['g'] * (p['betaS'] * (S[1:] - S[:-1]) - p['betaT'] * (T[1:] - T[:-1]))

    # Homogenize layers with negative buoyancy
    if np.any(B < 0):
        # Find the indices of negative buoyancy entries
        inx = np.where(B < 0)[0]  # Indices of layers with negative buoyancy
        for k in inx:
            # Indices of the unstable pair of layers
            inds = [k, k + 1]

            # Homogenize temperature and salinity (volume-weighted average)
            T[inds] = np.sum(T[inds] * V[inds]) / np.sum(V[inds])
            S[inds] = np.sum(S[inds] * V[inds]) / np.sum(V[inds])

    # Put homogenized solution back into the output structure
    s['T'][:, i] = T
    s['S'][:, i] = S

    return s


def step_solution_forwards(i, p, s):
    """
    Perform an Euler timestep for the simulation.

    Parameters:
    - i: Current timestep index (integer)
    - p: Dictionary containing model parameters
    - s: Dictionary containing state variables at timestep i (T, S, etc.)

    Returns:
    - s: Updated dictionary with state variables at timestep i+1
    """

    # Compute the tracer variables at timestep i+1
    if s['QTp'].ndim == 1:
        # If QTp is 1D, update temperature and salinity directly
        s['T'][:, i+1] = s['T'][:, i] + s['dt'][i] * p['sid'] * (s['QTp'][:, i] + s['QTs'][:, i] + s['QTk'][:, i] + s['QTi'][:, i] + s['QTv'][:, i]) / s['V']
        s['S'][:, i+1] = s['S'][:, i] + s['dt'][i] * p['sid'] * (s['QSp'][:, i] + s['QSs'][:, i] + s['QSk'][:, i] + s['QSi'][:, i] + s['QSv'][:, i]) / s['V']
    else:
        # If QTp is 3D (i.e., more than one plume), it will sum the fluxes across the plume dimension
        s['T'][:, i+1] = s['T'][:, i] + s['dt'][i] * p['sid'] * (np.sum(s['QTp'][:, :, i],axis=0) + s['QTs'][:, i] + s['QTk'][:, i] + s['QTi'][:, i] + s['QTv'][:, i]) / s['V']
        s['S'][:, i+1] = s['S'][:, i] + s['dt'][i] * p['sid'] * (np.sum(s['QSp'][:, :, i],axis=0) + s['QSs'][:, i] + s['QSk'][:, i] + s['QSi'][:, i] + s['QSv'][:, i]) / s['V']
    
    return s


def check_solution(i, s):
    """
    Check for errors in volume conservation at timestep i.

    Parameters:
    - i: Current timestep index (integer)
    - s: Dictionary containing state variables and fluxes at timestep i

    Returns:
    - status: Integer (0 if volume is conserved, 1 if there's a volume error)
    """

    status = 0  # Initialize status to 0 (no error)

    # Check that volume is conserved
    #if s['QVp'].ndim == 1:
    #    # If QVp is 1D (single plume), compute max dV/dt
    #    maxdVdt = np.max(np.abs(s['QVp'][:, :, i].sum(axis=0) + s['QVs'][:, i] + s['QVk'][:, i] + s['QVi'][:, i] + s['QVv'][:, i]))
    #else:
    #    # If QVp is 2D (multiple plumes), sum the fluxes across plumes
    maxdVdt = np.max(np.abs(np.sum(s['QVp'][:, :, i], axis=0) + s['QVs'][:, i] + s['QVk'][:, i] + s['QVi'][:, i] + s['QVv'][:, i]))

    # If the maximum change in volume is too large, flag an error
    if maxdVdt > 1e-8:
        print(f"Warning: volume possibly not conserved, max dV/dt = {maxdVdt}")
        status = 1  # Set status to 1 to indicate an error

    return status


def get_final_output(p, t, s, status):
    """
    Get the final output of the simulation, saving the necessary data at
    the specified time steps.

    Parameters:
    - p: Parameter dictionary
    - t: Time vector
    - s: Solution dictionary
    - status: Simulation status

    Returns:
    - s: Updated solution dictionary with final outputs
    """
    
    # Recompute fluxes for the final output
    s = ffl.compute_fluxes(s['T'].shape[1]-1, p, s)
    
    # If there's an error, ensure that all time steps are properly saved
    if status == 1:
        p['t_save'] = t

    # Get the indices of the timesteps to save
    inx = np.isin(t, p['t_save'])

    # Update solution with the selected time steps
    s['t'] = t[inx]
    s['T'] = s['T'][:, inx]
    s['S'] = s['S'][:, inx]

    # Calculate the vertical grid
    ints = np.concatenate(([0], np.cumsum(s['H'])))
    s['z'] = 0.5 * (ints[:-1] + ints[1:])

    # Plume exchange fluxes
    s['QVp'] = s['QVp'][:, :, inx]
    s['QTp'] = s['QTp'][:, :, inx]
    s['QSp'] = s['QSp'][:, :, inx]
    s['QMp'] = s['QMp'][:, :, inx]
    s['QEp'] = s['QEp'][:, :, inx]
    s['knb'] = s['knb'][:, inx]
    
    # Shelf exchange
    s['QVs'] = s['QVs'][:, inx]
    s['QTs'] = s['QTs'][:, inx]
    s['QSs'] = s['QSs'][:, inx]
    s['Ss'] = s['Ss'][:, inx]
    s['Ts'] = s['Ts'][:, inx]
    s['phi'] = s['phi'][:, inx]

    # Vertical mixing fluxes
    s['QVk'] = s['QVk'][:, inx]
    s['QTk'] = s['QTk'][:, inx]
    s['QSk'] = s['QSk'][:, inx]

    # Vertical fluxes
    s['QVv'] = s['QVv'][:, inx]
    s['QTv'] = s['QTv'][:, inx]
    s['QSv'] = s['QSv'][:, inx]

    # Iceberg fluxes 
    s['QVi'] = s['QVi'][:, inx]
    s['QTi'] = s['QTi'][:, inx]
    s['QSi'] = s['QSi'][:, inx]
    s['QMi'] = s['QMi'][:, inx]
    
    # melt rates
    s['plumemeltrate']   = np.zeros(s['QVp'].shape)
    s['icebergmeltrate'] = np.zeros(s['QMi'].shape)
    for k in range(len(s['t'])):
        for j in range(len(p['wp'])):
            s['plumemeltrate'][j, :, k] = p['sid'] * s['QMp'][j, :, k] / (p['wp'][j] * s['H'])

        with np.errstate(divide='ignore',invalid='ignore'): # we handle this possible division by zero down below
            iceberg_meltrate = p['sid'] * s['QMi'][:,k] / s['I']
        iceberg_meltrate[s['I'] == 0] = 0  # Set to 0 where no icebergs
        s['icebergmeltrate'][:,k] = iceberg_meltrate

    # Subglacial discharge
    s['Qsg'] = s['Qsg'][:, inx]

    # Store the final status
    s['status'] = status
    
    return s



def run_model(p, t, f, a):
    """
    RUN_MODEL Run fjordRPM simulation.
    s = RUN_MODEL(p, t, f, a, path_out) runs the fjordRPM simulation for
    parameters structure p, time vector t, forcings structure f, initial
    conditions structure a, and returns solution structure s.
    """

    # Check for errors in the given inputs
    status = check_inputs(p, t, f, a)

    # Preallocate and initialise variables
    s = initialise_variables(p, t, f, a)

    # The timestepping loop
    for i in range(len(t) - 1):
        # Homogenize heat & salt content of layers where density is unstable
        s = homogenise_unstable_layers(i, p, s)

        # Compute the fluxes ready for next timestep
        s = ffl.compute_fluxes(i, p, s)

        # Step the tracers forward in time
        s = step_solution_forwards(i, p, s)


        # Check for errors at this new timestep
        status = check_solution(i, s)

    # Subsample solution structure to requested output frequency
    s = get_final_output(p, t, s, status)

    return s
