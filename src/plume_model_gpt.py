#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:28:36 2024

@author: original code from Slater et al. (2025; GMD), translated to Python by ChatGPT and reviewed by Martim Mas e Braga
"""
import numpy as np

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
    
    # Locally allocate constants that don't change to speed up computation
    ci = p['ci']
    cw = p['cw']
    GS = p['GS']
    GT = p['GT']
    l1 = p['l1']
    l2 = p['l2']
    l3 = p['l3']
    Ti = p['Ti']
    l = p['l']
    Cd = p['Cd']
    
    # Precompute what is used repeatedly to speed up computation
    l2_plus_l3_z0 = l2 + l3 * z0 
    
    # Precalculate the terms for the quadratic equation to speed up computation
    a1 = l1 * (ci * GS - cw * GT)
    a2 = cw * GT * (T0 - l2 - l3 * z0) + GS * (ci * (l2_plus_l3_z0 - l1 * S0 - Ti) + l)
    a3 = -GS * S0 * (ci * (l2_plus_l3_z0 - Ti) + l)
    discriminant = a2**2 - 4 * a1 * a3 # Calculate the discriminant (inside the square root)
    
    Sb = (-a2 + np.sqrt(discriminant)) / (2 * a1) # Calculate the salinity at the boundary (Sb)
    Tb = l1 * Sb + l2_plus_l3_z0 # Calculate the temperature at the boundary (Tb)
    
    # Calculate the melt rate (mdot)
    mdot = cw * np.sqrt(Cd) * GT * u0 * (T0 - Tb) / (l + ci * (Tb - Ti))
    
    return mdot, Tb, Sb




def run_plume(j, p, kgl, H0, S0, T0, Qsg0):
    """
    RUN_PLUME computes the plume dynamics for a given plume index and other model parameters.
    """
    # Reverse the properties to orient deepest layers first
    H0 = np.flipud(H0)
    S0 = np.flipud(S0)
    T0 = np.flipud(T0)
    kgl = len(H0) - 1 - kgl  # Adjust grounding line index

    # Initialize variables
    n = len(H0)
    mdot, Tb, Sb   = np.zeros(n), np.zeros(n), np.zeros(n)
    Qent, Qmelt    = np.zeros(n), np.zeros(n)
    Tp, Sp, gp     = np.zeros(n), np.zeros(n), np.zeros(n)
    b, u, edot     = np.zeros(n), np.zeros(n), np.zeros(n)
    QV, QM, QT, QS = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)

    # Set properties at the grounding line (initial values)
    Tp[kgl] = p['l2'] + p['l3'] * p['Hgl'][j]
    Sp[kgl] = 0
    gp[kgl] = p['g'] * (p['betaS'] * (S0[kgl] - Sp[kgl]) - p['betaT'] * (T0[kgl] - Tp[kgl]))
    b[kgl] = (p['alphap'] * (Qsg0 / p['wp'][j])**2 / gp[kgl])**(1/3)
    u[kgl] = Qsg0 / (p['wp'][j] * b[kgl])
    edot[kgl] = p['alphap'] * u[kgl]
    QV[kgl] = u[kgl] * b[kgl]
    QM[kgl] = u[kgl]**2 * b[kgl]
    QT[kgl] = b[kgl] * u[kgl] * Tp[kgl]
    QS[kgl] = b[kgl] * u[kgl] * Sp[kgl]

    mdot[kgl], Tb[kgl], Sb[kgl] = meltrate(p, u[kgl], Tp[kgl], Sp[kgl], p['Hgl'][j])

    # Loop over layers to compute plume dynamics
    k = kgl
    while gp[k] > 0 and k < n - 1:
        k += 1
        dz = p['Hgl'][j] - np.sum(H0[:kgl + 1]) if k == kgl + 1 else H0[k - 1]
        
        # Compute fluxes for current layer
        Qent[k - 1] = dz * edot[k - 1]
        Qmelt[k - 1] = dz * mdot[k - 1]
        QV[k] = QV[k - 1] + Qent[k - 1] + Qmelt[k - 1]
        QM[k] = QM[k - 1] + dz * (b[k - 1] * gp[k - 1] - p['Cd'] * u[k - 1]**2)
        QT[k] = QT[k - 1] + dz * (edot[k - 1] * T0[k - 1] + mdot[k - 1] * Tb[k - 1] - p['Cd']**0.5 * p['GT'] * u[k - 1] * (Tp[k - 1] - Tb[k - 1]))
        QS[k] = QS[k - 1] + dz * (edot[k - 1] * S0[k - 1])

        # Update plume model variables
        b[k] = QV[k]**2 / QM[k]
        u[k] = QM[k] / QV[k]
        Tp[k] = QT[k] / QV[k]
        Sp[k] = QS[k] / QV[k]
        gp[k] = p['g'] * (p['betaS'] * (S0[k] - Sp[k]) - p['betaT'] * (T0[k] - Tp[k]))
        edot[k] = p['alphap'] * u[k]
        mdot[k], Tb[k], Sb[k] = meltrate(p, u[k], Tp[k], Sp[k], p['Hgl'][j] - np.sum(H0[kgl:k]))

    # Scale volume fluxes for plume width and reverse order to orient shallowest first
    Qent = np.flipud(Qent) * p['wp'][j]
    Qmelt = np.flipud(Qmelt) * p['wp'][j]

    # Find the neutral buoyancy index
    knb = len(H0) - len(gp[gp != 0])

    return Qent, Qmelt, knb



