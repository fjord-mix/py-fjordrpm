#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:28:36 2024

@author: original code from Slater et al. (2025; GMD), translated to Python by ChatGPT and reviewed by Martim Mas e Braga
"""
import numpy as np

def meltrate(p, u0, T0, S0, z0):
    """
    MELTRATE computes the submarine melt rate, temperature, and salinity.
    
    Parameters:
    - p: Dictionary containing model parameters (e.g., `l1`, `ci`, `cw`, etc.)
    - u0: Flow velocity (m/s) at the given depth
    - T0: Temperature at the given depth (°C or K depending on model units)
    - S0: Salinity at the given depth (psu or equivalent)
    - z0: Depth at which melt is calculated (m)

    Returns:
    - mdot: Melt rate (m/s)
    - Tb: Temperature at the boundary (°C or K)
    - Sb: Salinity at the boundary (psu or equivalent)
    """
    
    # Calculate the coefficients for the quadratic equation
    a1 = p['l1'] * (p['ci'] * p['GS'] - p['cw'] * p['GT'])
    a2 = p['cw'] * p['GT'] * (T0 - p['l2'] - p['l3'] * z0) + p['GS'] * (p['ci'] * (p['l2'] + p['l3'] * z0 - p['l1'] * S0 - p['Ti']) + p['l'])
    a3 = -p['GS'] * S0 * (p['ci'] * (p['l2'] + p['l3'] * z0 - p['Ti']) + p['l'])
    
    # Calculate the salinity at the boundary
    Sb = (-a2 + np.sqrt(a2**2 - 4 * a1 * a3)) / (2 * a1)
    
    # Calculate the temperature at the boundary
    Tb = p['l1'] * Sb + p['l2'] + p['l3'] * z0
    
    # Calculate the melt rate
    mdot = p['cw'] * p['Cd']**(1/2) * p['GT'] * u0 * (T0 - Tb) / (p['l'] + p['ci'] * (Tb - p['Ti']))
    
    return mdot, Tb, Sb


def run_plume(j, p, kgl, H0, S0, T0, Qsg0):
    """
    RUN_PLUME computes the plume dynamics for a given plume index and other model parameters.

    Parameters:
    - j: The plume index
    - p: Dictionary of model parameters (e.g., physical constants, plume properties)
    - kgl: Index of the grounding line in the depth profile
    - H0: Array of initial layer thicknesses
    - S0: Array of initial salinities at the model layers
    - T0: Array of initial temperatures at the model layers
    - Qsg0: Subglacial discharge for this plume at the current timestep

    Returns:
    - Qent: Entrainment flux for each layer
    - Qmelt: Melt flux for each layer
    - knb: Index of the neutral buoyancy layer
    """
    # Reverse the properties to orient deepest layers first
    Ta = np.flipud(T0)
    Sa = np.flipud(S0)
    H0 = np.flipud(H0)
    kgl = len(H0)-1 - kgl  # Adjust grounding line index
    mdot,Tb,Sb = [np.zeros(len(H0)) for _ in range(3)]

    # Initialise output variables
    Qent = np.zeros(len(H0))
    Qmelt = np.zeros(len(H0))

    # Initial plume model variables
    Tp = np.zeros(len(H0))
    Sp = np.zeros(len(H0))
    gp = np.zeros(len(H0))
    b = np.zeros(len(H0)) 
    u = np.zeros(len(H0)) 
    edot = np.zeros(len(H0))
    QV = np.zeros(len(H0))
    QM = np.zeros(len(H0))
    QT = np.zeros(len(H0))
    QS = np.zeros(len(H0))

    # Set properties at the grounding line (initial values)
    Tp[kgl] = p['l2'] + p['l3'] * p['Hgl'][j]
    Sp[kgl] = 0
    gp[kgl] = p['g'] * (p['betaS'] * (Sa[kgl] - Sp[kgl]) - p['betaT'] * (Ta[kgl] - Tp[kgl]))
    b[kgl] = (p['alphap'] * (Qsg0 / p['wp'][j])**2 / gp[kgl])**(1/3)
    u[kgl] = Qsg0 / (p['wp'][j] * b[kgl])
    edot[kgl] = p['alphap'] * u[kgl]
    QV[kgl] = u[kgl] * b[kgl]
    QM[kgl] = u[kgl]**2 * b[kgl]
    QT[kgl] = b[kgl] * u[kgl] * Tp[kgl]
    QS[kgl] = b[kgl] * u[kgl] * Sp[kgl]
    mdot[kgl], Tb[kgl], Sb[kgl] = meltrate(p, u[kgl], Tp[kgl], Sp[kgl], p['Hgl'][j])

    # Loop over layers to compute plume dynamics
    # TODO: can we make this more efficient?
    k = kgl
    while gp[k] > 0 and k < len(H0)-1:
        # Advance the fluxes
        k += 1
        if k == kgl + 1:
            dz = p['Hgl'][j] - (p['H'] - np.sum(H0[:kgl+1]))
        else:
            dz = H0[k - 1]

        Qent[k - 1] = dz * edot[k - 1]
        Qmelt[k - 1] = dz * mdot[k - 1]
        QV[k] = QV[k - 1] + Qent[k - 1] + Qmelt[k - 1]
        QM[k] = QM[k - 1] + dz * (b[k - 1] * gp[k - 1] - p['Cd'] * u[k - 1]**2)
        QT[k] = QT[k - 1] + dz * (edot[k - 1] * Ta[k - 1] + mdot[k - 1] * Tb[k - 1] - p['Cd']**0.5 * p['GT'] * u[k - 1] * (Tp[k - 1] - Tb[k - 1]))
        QS[k] = QS[k - 1] + dz * (edot[k - 1] * Sa[k - 1])

        # Update plume model variables
        b[k] = QV[k]**2 / QM[k]
        u[k] = QM[k] / QV[k]
        Tp[k] = QT[k] / QV[k]
        Sp[k] = QS[k] / QV[k]
        gp[k] = p['g'] * (p['betaS'] * (Sa[k] - Sp[k]) - p['betaT'] * (Ta[k] - Tp[k]))
        edot[k] = p['alphap'] * u[k]
        mdot[k], Tb[k], Sb[k] = meltrate(p, u[k], Tp[k], Sp[k], p['Hgl'][j] - np.sum(H0[kgl:k]))

    # Scale volume fluxes for plume width and reverse order to orient shallowest first
    Qent = np.flipud(Qent) * p['wp'][j]
    Qmelt = np.flipud(Qmelt) * p['wp'][j]

    # Find the neutral buoyancy index
    knb = len(H0) - len(gp[gp != 0]) #

    return Qent, Qmelt, knb


