#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 7 21:54:56 2020
@author: Nicholas Vieira
@profiles.py

Compute density and temperature profiles applicable to a kilonova on a radial 
grid. Profiles are taken from:
--> Bulla 19, MNRAS 489, 5037-5045 (their equations 2, 3, 12)

"""

import numpy as np

#### DENSITY PROFILES #########################################################
def rho_init(r, a=1, beta=-3):
    """
    r: distance from center in cm
    a: scaling factor to get a desired ejecta mass (optional; default 1)
    beta: slope of density power-law (optional; default -3)
    
    Compute initial density profile assuming some scaling factor a and a power
    law slope of -3.
    """
    return a * r**beta

def rho_init_Kasen13(t, beta_ej=0.2, Mej=0.04):
    """
    t: time since homology set in the ejecta (effectively time post-merger 
    for a KNe) in ***days***
    beta_ej: characteristic ejecta velocity v/c (optional; default 0.2)
    Mej: ejecta mass in solar masses (optional; default 0.04)
    
    Estimate the density at ~1 day post-merger for a BNS KNe, as given in: 
    --> Kasen+13, ApJ 774, 1 (their equation 3)
    
    Note that this density is just the density at the photosphere at this time.
    """
    return (2.8e-13) * Mej/0.01 * (beta_ej/0.1)**-3 * t**-3


def rho(r, t, t0=1.5, beta=-3, rho_i=None, a=1):
    """
    r: distance from center in cm
    t: time post-merger in days
    t0: time at which rho_init was computed in days (optional; default 1.5)
    beta: slope of density power-law (optional; default -3)
    rho_i: initial density at <t0> days (optional; default None, which uses 
           function rho_init() above to estimate)
    a: scale factor for computing rho_init() (optional; default 1; only 
       relevant if rho_i is not given)
    
    Compute time and space-dependent density profile.
    """
    # if no initial density assigned
    if type(rho_i) == type(None):
        if type(t) in [list, np.ndarray]: # if a list of times is given
            ret = [rho_init(r, a=a) * (tim/t0)**beta for tim in t]
            return np.array(ret)
        else: # if a single time is given
            return rho_init(r, a=a) * (t/t0)**beta
        
    # if an initial density is assigned
    else:
        if type(t) in [list, np.ndarray]: # if a list of times is given            
            ret = [rho_i * (tim/t0)**beta for tim in t]
            return np.array(ret)
        else: # if a single time is given
            return rho_i * (t/t0)**beta


#### TEMPERATURE PROFILES #####################################################
def T(t, t0=1.5, alpha=-0.4, T_init=5000, T_floor=0):
    """
    t: time post-merger in days
    t0: time at which T_init was computed in days (optional; default 1.5)
    alpha: slope of temperature power-law (optional; default -0.4)
    T_init: temperature at <t0> days (optional; default 5000K)
    T_floor: temperature minimum (optional; default 0) 
    
    Compute time-dependent temperature profile. No spatial dependence.
    """
    
    if type(t) in (list, np.ndarray): # if a list of times is given
        T = [max(T_init * (tim/t0)**alpha, T_floor) for tim in t]
        return np.array(T)
    else: # if a single time is given
        return max([T_init * (t/t0)**alpha, T_floor])
