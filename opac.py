#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 7 17:12:42 2020
@author: Nicholas Vieira
@opac.py

Compute the electron scattering opacity, or compute the wavelength-dependent 
bound-bound opacities for either a lanthanide-rich ejecta or lanthanide-poor 
ejecta. 

Bound-bound opacities are taken from:
--> Bulla 19, MNRAS 489, 5037-5045 (their figure 2)

which in turn used the following paper for their calculations:
--> Tanaka+19, arXiv:1906.08914
"""

TREF = 1.5 # reference time, in days
GAMMA = 1 # gamma power-law parameter

#### ELECTRON SCATTERING OPACITIES ############################################
def kappa_es(t):
    """
    t: time in days
    """
    return 0.01 * (t/TREF)**(-GAMMA)


#### LANTHANIDE-POOR BOUND-BOUND OPACITIES ####################################
def prefac_bb_lanthpoor(wavelen):
    """
    wavelen: photon wavelength in microns
    
    lanthanide-poor bound-bound opacities are described by:
    log10(k_bb) = 2 - wavelen*4.5;  wavelen < 1 micron
    k_bb = 5e-3;                    wavelen > 1 micron
    """    
    if wavelen <= 1:
        return 10**(2 - wavelen*4.3)
    else:
        return 5e-3


def kappa_bb_lanthpoor(t, prefac=5e-3):
    """
    t: time in days
    prefac: opacity at t=1.5 days for some wavelength
    """
    return prefac * (t/TREF)**GAMMA


#### LANTHANIDE-RICH BOUND-BOUND OPACITIES #################################### 
def prefac_bb_lanthrich(wavelen):
    """
    wavelen: photon wavelength in microns
    
    lanthanide-rich bound-bound opacities are described by:
    log(k_bb) = 2 - wavelen*2; wavelen < 1 micron
    k_bb = 1;                  wavelen > 1 micron
    """    
    if wavelen <= 1:
        return 10**(2 - wavelen*2)
    else:
        return 1


def kappa_bb_lanthrich(t, prefac=1.0):
    """
    t: time in days
    prefac: opacity at t=1.5 days for some wavelength
    """
    return prefac * (t/TREF)**GAMMA
