#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 7 17:13:24 2020
@author: Nicholas Vieira
@distribs.py

Generate a Maxwell-Boltzmann electron distribution, a Planck photon 
distribution, or a Planck function multiplied by the total frequency-dependent
opacity of the system.
"""

import numpy as np
from astropy.constants import k_B, c, m_e, h
import opac

k_B = k_B.cgs.value # convert to cgs, take value
c = c.cgs.value    
m_e = m_e.cgs.value
h = h.cgs.value

def gen_maxwell_boltzmann(T, n=10000):
    """
    T: electron distribution temperature in K
    n: number of samples (optional; default 10000)
    
    Compute the Maxwell-Boltzmann distribution for some T. Then, take <n>
    random samples from this distribution.

    The distribution from which the function draws will have max(n, 10000)
    values from which to draw. 
    
    Returns an array of random beta=v/c drawn from the Maxwell-Boltzmann 
    distribution. 
    """
    # compute distribution
    prefac = (m_e*0.5/(np.pi*k_B*T))**1.5 * 4 * np.pi
    v = np.linspace(0, c, max(10000, n)) # compute for beta = 0 - 1
    dist = prefac * v**2 * np.exp(-m_e*v**2/(2*k_B*T))
    dist = dist/np.sum(dist) # normalize 
    
    # random draws weighted by probability distribution
    draws = np.random.choice(np.linspace(0, 1, max(10000, n)), size=n, p=dist)
    
    return draws


def gen_Planck(T, n=10000):
    """
    T: photon distribution temperature in K
    n: number of samples (optional; default 10000)
    
    Compute the Planck distribution for some T. Then, take <n> random samples 
    from this distribution.

    The distribution from which the function draws will have max(n, 10000)
    values from which to draw. 
   
    Returns an array of random photon **energies** drawn from the Planck 
    distribution, in cgs units.
    """    
    # use Wien's displacement law to find suitable range of frequencies for kT
    nu_max = 2.82*k_B*T/h
    # compute for nu = nu_max/10000 to nu_max*10000
    nu = np.logspace(np.log10(nu_max)-4, np.log10(nu_max)+4, max(10000, n)) 
    
    # finally, compute the distribution
    dist = (2 * h * nu**3 / c**2)/(np.expm1(h*nu/(k_B*T)))
    dist = dist/np.sum(dist) # normalize 
    
    # random draws weighted by probability distribution
    draws = np.random.choice(np.logspace(np.log10(nu_max)-4, 
                                         np.log10(nu_max)+4, max(10000, n)), 
                             size=n, p=dist)
    
    return draws*h


def gen_thermal_emissivity(T, t, Ye, n=10000):
    """
    T: photon distribution temperature in K
    t: time post-merger in days
    Ye: electron fraction (< 0.25 = lanthanide-rich, > 0.25 = lanthanide-poor)
    n: number of samples (optional; default 10000)
    
    Compute the Planck distribution for some T, **multiplied** by kappa_tot, 
    which is frequency-dependent and represents the total opacity of the 
    system. Then, take <n> random samples from this distribution.
    
    The distribution from which the function draws will have max(n, 10000)
    values from which to draw. 
    
    Returns an array of random photon **energies** drawn from this thermal 
    emissivity distribution, in cgs units.
    """    
    # use Wien's displacement law to find suitable range of frequencies for kT
    nu_max = 2.82*k_B*T/h
    # compute for nu = nu_max/10000 to nu_max*10000
    nu = np.logspace(np.log10(nu_max)-4, np.log10(nu_max)+4, max(10000, n)) 
    wavelens = c/nu * 1e4 # wavelength in microns    
    
    # finally, compute the distribution
    if Ye < 0.25:
        prefac = np.array([opac.prefac_bb_lanthrich(w) for w in wavelens])
        total_opac = opac.kappa_es(t) + opac.kappa_bb_lanthrich(t, prefac)
    else:
        prefac = np.array([opac.prefac_bb_lanthpoor(w) for w in wavelens])
        total_opac = opac.kappa_es(t) + opac.kappa_bb_lanthpoor(t, prefac)
        
    dist = (2 * h * nu**3 / c**2)/(np.expm1(h*nu/(k_B*T))) * total_opac
    dist = dist/np.sum(dist) # normalize 
    
    # random draws weighted by probability distribution
    draws = np.random.choice(np.logspace(np.log10(nu_max)-4, 
                                         np.log10(nu_max)+4, max(10000, n)), 
                             size=n, p=dist)
    
    return draws*h