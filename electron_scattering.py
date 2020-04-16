#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 22:01:37 2020
@author: Nicholas Vieira
@electron_scattering.py

Compute outgoing angles and energies for a photon scattering a single time off 
of an electron.
"""

import numpy as np

def y_gen_uniform(beta, n=10000):
    """
    beta: ejecta v/c
    n: no. of samples (optional; default 10,000)
    
    generate uniform distribution of y with <n> samples between 
    (1 - beta**2)/(2*beta) and (1 + beta**2)/(2*beta) 
    applying the 'transformation method'
    """
    return np.linspace(((1 - beta)**2)/(2*beta), # min
                       ((1 + beta)**2)/(2*beta), # max
                       n)
    
def cosalpha_gen(n=10000):
    """
    n: no. of samples (optional; default 10,000)
    
    generate correct distribution for cos(alpha) with <n> samples
    applying the 'rejection method'
    """

    # uniform distributions for 2 random variates
    x_uni = np.linspace(0, 2, n)
    y_uni = np.linspace(-1, 1, n)
    
    # perform rejection method until cosalpha distribution has <n> values 
    cosalpha_dist = []
    while len(cosalpha_dist) < n:
        x_rand = np.random.choice(x_uni) # pick random x
        y_rand = np.random.choice(y_uni) # pick random y
        if x_rand < 1 + y_rand**2: # if satisfied, add y to cosalpha distrib
            cosalpha_dist.append(y_rand)
        # if not satisfied, reject and try again
        
    return cosalpha_dist


def scatter(beta, mu_in, E_in):
    """
    beta: ejecta velocity v/c
    mu_in: ingoing photon propagation angle wrt to the velocity vector of the 
           ejecta
    E_in: rest-frame energy of the incident photon
    
    For a photon packet with energy E_in travelling with initial propagation 
    angle cosine mu_in which collides with an electron with some beta = v/c,
    compute the outgoing propagation angle cosine mu_out and the outgoing 
    energy E_out after the scattering.
    """
    # sample cosalpha using rejection method
    cosalpha = np.random.choice(cosalpha_gen(1)) # just one
    alpha = np.arccos(cosalpha)        
    # sample phi uniformly from 0 to 2pi
    phi = np.random.choice(np.linspace(0, 2*np.pi, 1000))
    
    # incident angle in co-moving frame
    mu_in_cmf = (mu_in - beta)/(1 - beta*mu_in)
    theta_in_cmf = np.arccos(mu_in_cmf) 
    
    # outgoing angle in co-moving frame
    mu_out_cmf = cosalpha*mu_in_cmf 
    mu_out_cmf -= np.sin(alpha)*np.sin(theta_in_cmf)*np.cos(phi)
    
    # outgoing angle in rest frame
    mu_out = (mu_out_cmf + beta)/(1 + mu_out_cmf*beta)
    
    # compute the outgoing packet energy
    E_out = E_in * (1 - mu_in*beta) / (1 - mu_out*beta)
    
    return mu_out, E_out
    
    
            
        
    


