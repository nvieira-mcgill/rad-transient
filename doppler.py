#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 17:18:50 2020
@author: Nicholas Vieira
@doppler.py

Computing the Lorentz factor from some beta = v/c. Doppler shifts for 
transforming wavelength from rest-frame to co-moving frame.
"""

def gamma(beta):
    return (1 - beta**2)**(-0.5)

def wavelen_cmf(wavelen, beta, mu):
    """
    wavelen: photon rest-frame wavelength
    beta: ejecta v/c
    mu: cosine of incident angle of photon 
    """
    return wavelen * (1 - mu*beta)**(-1)

