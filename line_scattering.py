#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 23:13:04 2020
@author: Nicholas Vieira
@line_scattering.py

Scatter a photon packet on a known atomic line. A crude model for an atomic 
line which invokes the Sobolev approximation. Based on the work of:
--> Kasen et al. 2006, ApJ 651, 366-380

"""

import numpy as np

def scatter(tau_l, beta_ej, mu_in, E_in):
    """
    tau_l: optical depth of the line
    beta_ej: velocity v/c of the ejecta in which the interaction occurs
    mu_in: ingoing photon propagation angle wrt to the velocity vector of the 
           ejecta
    E_in: rest-frame energy of the incident photon
    
    Using the optical depth of the line, determine whether the photon packet is
    trapped in the line. If trapped, scatter the photon in the line until it 
    escapes.
    """
    
    # compute interaction probability, compare to random draw in [0, 1)
    z = np.random.choice(np.linspace(0, 0.9999, 1000))
    if z < 1 - np.exp(-tau_l): # an interaction occurs
        
        # compute escape probability 
        p_escape = (1 - np.exp(-tau_l)) / tau_l
        
        x = 0 # keep scattering until photon escapes resonance
        while x < p_escape: 
            print("line scatter!")
            # randomly select outgoing scattering angle (isotropic)
            mu_out = np.random.choice(np.linspace(0, 1.0, 1000))
            # compute the outgoing packet energy
            E_out = E_in * (1 - mu_in*beta_ej) / (1 - mu_out*beta_ej)        
            
            # generate new x, see if photon escapes
            x = np.random.choice(np.linspace(0, 1.0, 1000))
            
        return mu_out, E_out 
    
    else: # photon does not fall into resonance with line 
        return mu_in, E_in
