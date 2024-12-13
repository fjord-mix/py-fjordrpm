#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:26:32 2024

@author: original code from Slater et al. (2025; GMD), translated to Python by ChatGPT and reviewed by Martim Mas e Braga
"""
import numpy as np
import plume_model as pl
import warnings

# we use that to supress division by zero warnings in get_mixing_fluxes and get_iceberg_fluxes
def fxn():
    warnings.warn("Runtime", RuntimeWarning)


def get_plume_fluxes(i, p, s):
    """
    GET_PLUME_FLUXES computes the plume fluxes for the given parameters and solution at timestep i.

    Parameters:
    - i: Current timestep index
    - p: Dictionary containing model parameters (e.g., physical constants, plume properties)
    - s: Dictionary containing the solution structure (e.g., tracer concentrations, fluxes)

    Returns:
    - QVp0, QTp0, QSp0: Plume volume, temperature, and salinity fluxes for each plume at timestep i
    - QEp0, QMp0: Entrainment and melt fluxes for each plume at timestep i
    - knb0: Neutral buoyancy layer index for each plume at timestep i
    """
    # Get tracer variables at timestep i
    H0 = s['H']
    T0 = s['T'][:, i]
    S0 = s['S'][:, i]

    # Initialise outputs
    QVp0, QTp0, QSp0, QEp0, QMp0 = [np.zeros((len(p['wp']), p['N'])) for _ in range(5)]
    knb0 = np.zeros(len(p['wp']))

    # Loop over number of plumes
    for j in range(len(p['wp'])):
        # Get subglacial discharge at timestep i
        Qsg0 = s['Qsg'][j, i]
        kgl = s['kgl'][j]

        if Qsg0 != 0:  # If there is a plume
            # Plume dynamics
            if (i - 1) % p['run_plume_every'] == 0 or s['Qsg'][j, i - 1] == 0:  # If a plume update timestep
                Qent, Qmelt, knb = pl.run_plume(j, p, kgl, H0, S0, T0, Qsg0)
            else:  # Otherwise, use dynamics from previous time step
                Qent = s['QEp'][j, :, i - 1]
                Qmelt = s['QMp'][j, :, i - 1]
                knb = s['knb'][j, i - 1]

            # Compute fluxes in layers from grounding line to neutral buoyancy
            QVp0[j, knb:kgl] = -Qent[knb:kgl]
            QTp0[j, knb:kgl] = QVp0[j, knb:kgl] * T0[knb:kgl]
            QSp0[j, knb:kgl] = QVp0[j, knb:kgl] * S0[knb:kgl]

            # Compute fluxes into the neutral buoyancy layer
            Tsg0 = p['l2'] + p['l3'] * p['Hgl'][j]
            Teff = -p['l'] / p['cw']
            QVp0[j, knb] = Qsg0 + np.sum(Qmelt[knb:kgl]) - np.sum(QVp0[j, knb:kgl])
            QTp0[j, knb] = Qsg0 * Tsg0 + np.sum(Qmelt[knb:kgl]) * Teff - np.sum(QTp0[j, knb:kgl])
            QSp0[j, knb] = -np.sum(QSp0[j, knb:kgl])

            # Store entrainment, submarine melt flux, and neutral buoyancy
            QEp0[j, :] = Qent
            QMp0[j, :] = Qmelt
            knb0[j] = knb

    return QVp0, QTp0, QSp0, QEp0, QMp0, knb0


def get_shelf_fluxes(i, p, s):
    """
    GET_SHELF_FLUXES Compute shelf fluxes for the z-model.
    
    Parameters:
    - i: Current timestep (integer)
    - p: Dictionary containing model parameters (e.g., g, betaS, betaT, C0, etc.)
    - s: Dictionary containing solution structure at timestep i (e.g., H, T, S, Ts, Ss, etc.)

    Returns:
    - QVs0: Shelf volume fluxes (numpy array)
    - QTs0: Shelf heat fluxes (numpy array)
    - QSs0: Shelf salt fluxes (numpy array)
    - phi0: Potential energy across layers (numpy array)
    """
    
    # Get tracer variables at timestep i
    H0 = s['H']
    T0 = s['T'][:, i]
    S0 = s['S'][:, i]
    Ts0 = s['Ts'][:, i]
    Ss0 = s['Ss'][:, i]
    
    # Initialise variables
    phi0 = np.zeros(p['N'])
    QVs0 = np.zeros(p['N'])

    # Compute the fjord-to-shelf reduced gravity (gp)
    gp = p['g'] * (p['betaS'] * (Ss0[:s['ksill']] - S0[:s['ksill']]) -
                   p['betaT'] * (Ts0[:s['ksill']] - T0[:s['ksill']]))

    # Calculate the potentials over above-sill layers
    phi0[0] = gp[0] * H0[0] / 2
    for k in range(1, s['ksill']):
        phi0[k] = phi0[k - 1] + gp[k - 1] * H0[k - 1] / 2 + gp[k] * H0[k] / 2

    # Compute the above-sill fluxes before barotropic compensation (zero if C0=0)
    Q = p['C0'] * p['W'] * H0[:s['ksill']] * phi0[:s['ksill']] / p['L']

    # Compute the above-sill fluxes after barotropic compensation
    total_flux = np.sum(s['Qsg'][:, i]) + np.sum(s['QMp'][:, :, i])
    QVs0[:s['ksill']] = Q - H0[:s['ksill']] * (total_flux + np.sum(Q)) / np.sum(H0[:s['ksill']])

    # Compute the resulting heat and salt fluxes
    QTs0 = np.where(QVs0 >= 0, QVs0 * Ts0, QVs0 * T0)
    QSs0 = np.where(QVs0 >= 0, QVs0 * Ss0, QVs0 * S0)

    return QVs0, QTs0, QSs0, phi0


def get_mixing_fluxes(i, p, s):
    """
    Compute vertical tracer mixing fluxes at timestep i.

    Parameters:
    - i: Current timestep (integer)
    - p: Dictionary of model parameters (g, betaS, betaT, Kb, K0, Ri0, etc.)
    - s: Dictionary containing model state variables at timestep i (H, T, S, QVp, QVs, etc.)

    Returns:
    - QVk0: Vertical volume fluxes (numpy array)
    - QTk0: Vertical heat fluxes (numpy array)
    - QSk0: Vertical salt fluxes (numpy array)
    """
    with warnings.catch_warnings(action="ignore"):
        fxn()
        
    # Extract necessary variables at timestep i
    H0 = s['H']
    T0 = s['T'][:, i]
    S0 = s['S'][:, i]
    
    # Initialize vertical volume fluxes to zero (net mixing fluxes are zero)
    QVk0 = np.zeros_like(H0)
    
    # Compute the reduced gravity between adjacent layers
    gp = p['g'] * (p['betaS'] * (S0[1:] - S0[:-1]) - p['betaT'] * (T0[1:] - T0[:-1]))
    
    # Calculate horizontal velocity (u) from the volume fluxes
    if s['QVp'].ndim > 2:
        u = (np.sum(s['QVp'][:, :, i], axis=0) - s['QVs'][:, i]) / (2 * p['W'] * H0)
    else:
        u = (s['QVp'][:, i] - s['QVs'][:, i]) / (2 * p['W'] * H0)
    
    # Compute the Richardson number (Ri)
    du = u[1:] - u[:-1]
    Ri = gp * (H0[1:] + H0[:-1]) / (2 * du**2)
    Ri[du == 0] = p['Ri0']
    Ri[Ri > p['Ri0']] = p['Ri0']
    Ri[Ri < 0] = 0

    # Compute the diffusivity Kz as a function of the Richardson number
    Kz = p['Kb'] + ((Ri < p['Ri0']) & (Ri > 0)) * p['K0'] * (1 - (Ri / p['Ri0'])**2)**3
    
    # Compute the mixing fluxes for salinity (QS) and temperature (QT)
    QS = 2 * p['W'] * p['L'] * Kz * (S0[1:] - S0[:-1]) / (H0[1:] + H0[:-1])
    QT = 2 * p['W'] * p['L'] * Kz * (T0[1:] - T0[:-1]) / (H0[1:] + H0[:-1])
    
    # The final layer fluxes (QSk0, QTk0) are the net of the interface fluxes
    QTk0 = np.concatenate((QT, [0]),axis=0) - np.concatenate(([0], QT),axis=0)
    QSk0 = np.concatenate((QS, [0]),axis=0) - np.concatenate(([0], QS),axis=0)

    return QVk0, QTk0, QSk0


def get_iceberg_fluxes(i, p, s):
    """
    Compute iceberg fluxes at timestep i.

    Parameters:
    - i: Current timestep (integer)
    - p: Dictionary of model parameters (M0, l1, l2, l3, U0, etc.)
    - s: Dictionary containing model state variables at timestep i (H, T, S, I)

    Returns:
    - QVi0: Volume flux from icebergs (numpy array)
    - QTi0: Heat flux from icebergs (numpy array)
    - QSi0: Salt flux from icebergs (numpy array)
    - QMi0: Meltwater flux from icebergs (numpy array)
    """
    with warnings.catch_warnings(action="ignore"):
        fxn()
    
    # Extract necessary variables at timestep i
    H0 = s['H']
    T0 = s['T'][:, i]
    S0 = s['S'][:, i]
    I0 = s['I']
    
    # If there are no icebergs, return zero fluxes
    if p['M0'] == 0:
        QVi0, QTi0, QSi0, QMi0 = np.zeros_like(H0), np.zeros_like(H0), np.zeros_like(H0), np.zeros_like(H0)
    else:
        # Compute mean depth of boxes (zj)
        zj = np.cumsum(H0) - H0 / 2
        # Compute local freezing point Tf based on salinity and depth
        Tf = p['l1'] * S0 + p['l2'] + p['l3'] * zj
        # Compute melt flux into each box
        QMi0 = np.maximum(0, p['M0'] * (T0 - Tf) * I0)

        # Compute velocity scale and entrainment into upwelling
        Teff = -p['l'] / p['cw']  # Effective temperature of meltwater
        gmelt = p['g'] * (p['betaS'] * S0 - p['betaT'] * (T0 - Teff))  # Buoyancy difference
        vel = (QMi0 * gmelt * H0 / (p['alphai'] * I0))**(1/3)
        vel[(I0 == 0) | (gmelt <= 0)] = 0
        Qent = p['alphai'] * vel * I0

        # Compute length scale and fraction for upwelling
        gk = np.maximum(0, p['g'] * (p['betaS'] * (S0[1:] - S0[:-1]) - p['betaT'] * (T0[1:] - T0[:-1])))
        lice = (vel[1:]**2 / H0[1:]) * (H0[:-1] + H0[1:]) / gk
        lice[vel[1:] == 0] = 0
        fice = np.concatenate(([0], np.minimum(1, lice / H0[1:])))

        # Scale the upwelling flux
        Qentscaled = p['U0'] * fice * Qent

        # Compute the net volume, heat, and salt fluxes from upwelling
        QVi0 = np.concatenate((Qentscaled[1:], [0]),axis=0) - np.concatenate(([0], Qentscaled[1:]),axis=0)
        QTi0 = np.concatenate((Qentscaled[1:] * T0[1:], [0]),axis=0) - np.concatenate(([0], Qentscaled[1:] * T0[1:]),axis=0)
        QSi0 = np.concatenate((Qentscaled[1:] * S0[1:], [0]),axis=0) - np.concatenate(([0], Qentscaled[1:] * S0[1:]),axis=0)

        # Add the meltwater contributions to the heat and salt fluxes
        meltcont = np.concatenate((fice[1:], [0]),axis=0) * np.concatenate((QMi0[1:], [0]),axis=0) + (1 - fice) * QMi0
        QTi0 = QTi0 - meltcont * p['l'] / p['cw']
        QSi0 = QSi0 - meltcont * S0

    return QVi0, QTi0, QSi0, QMi0


def get_vertical_fluxes(i, s):
    """
    Compute vertical fluxes at timestep i.
    
    Parameters:
    - i: Current timestep (integer)
    - s: Dictionary containing model state variables at timestep i (T, S, QVp, QVs, QVi)
    
    Returns:
    - QVv0: Vertical volume fluxes (numpy array)
    - QTv0: Vertical heat fluxes (numpy array)
    - QSv0: Vertical salt fluxes (numpy array)
    """
    
    # Extract temperature and salinity at timestep i
    T0 = s['T'][:, i]
    S0 = s['S'][:, i]
    
    # Compute net flux imbalance (sum of volume, shelf, and iceberg fluxes)
    if s['QVp'].ndim == 2:  # If QVp has one plume
        Qnet = s['QVp'][:,i] + s['QVs'][:, i] + s['QVi'][:, i]
    else:  # If QVp has several plumes
        Qnet = np.sum(s['QVp'][:, :, i], axis=0) + s['QVs'][:, i] + s['QVi'][:, i]

    # Vertical flux required for no net volume change (cumulative sum)
    QVint = -np.cumsum(Qnet[:-1])

    # Compute temperature and salinity fluxes based on the direction of flux
    QTint = np.where(QVint < 0, QVint * T0[:-1], QVint * T0[1:])
    QSint = np.where(QVint < 0, QVint * S0[:-1], QVint * S0[1:])
    
    # Final vertical fluxes (flux going in and out of each layer)
    QVv0 = np.concatenate((QVint, [0]),axis=0) - np.concatenate(([0], QVint),axis=0)
    QTv0 = np.concatenate((QTint, [0]),axis=0) - np.concatenate(([0], QTint),axis=0)
    QSv0 = np.concatenate((QSint, [0]),axis=0) - np.concatenate(([0], QSint),axis=0)
    
    return QVv0, QTv0, QSv0


def compute_fluxes(i, p, s):
    """
    COMPUTE_FLUXES computes fluxes in the simulation at timestep i.
    This function calculates the plume fluxes, shelf fluxes, mixing fluxes,
    iceberg fluxes, and vertical fluxes for the given timestep.

    Parameters:
    - i: Current timestep index
    - p: Dictionary containing model parameters (e.g., physical constants, etc.)
    - s: Dictionary containing the solution state (e.g., tracer concentrations, fluxes, etc.)

    Returns:
    - s: Updated solution dictionary with computed fluxes for the current timestep
    """
    # Calculate plume fluxes
    s['QVp'][:,:,i], s['QTp'][:,:,i], s['QSp'][:,:,i], s['QEp'][:,:,i], s['QMp'][:,:,i], s['knb'][:,i] = get_plume_fluxes(i, p, s)

    # Calculate shelf fluxes
    s['QVs'][:,i], s['QTs'][:,i], s['QSs'][:,i], s['phi'][:,i] = get_shelf_fluxes(i, p, s)

    # Calculate tracer vertical mixing fluxes
    s['QVk'][:,i], s['QTk'][:,i], s['QSk'][:,i] = get_mixing_fluxes(i, p, s)

    # Calculate iceberg fluxes
    s['QVi'][:,i], s['QTi'][:,i], s['QSi'][:,i], s['QMi'][:,i] = get_iceberg_fluxes(i, p, s)

    # Calculate vertical fluxes
    s['QVv'][:,i], s['QTv'][:,i], s['QSv'][:,i] = get_vertical_fluxes(i, s)

    return s