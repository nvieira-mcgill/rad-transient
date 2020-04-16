#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 7 18:22:53 2020
@author: Nicholas Vieira
@heating.py

Compute radioactive heating rates, thermalization efficiencies, bolometric 
luminosities, and photospheric radii for a kilonova over time for known 
characteristic ejecta velocity, mass, and grey opacity.
"""

import numpy as np
from scipy.interpolate import interp2d
from scipy.integrate import trapz
from astropy.constants import M_sun, c, sigma_sb
Msun = M_sun.cgs.value
c = c.cgs.value
sigma_sb = sigma_sb.cgs.value

## constants in radioactive heating rate Lin
SIGMA_L_in = 0.11 # s
T0_L_in = 1.3 # s

def Lin(t, Mej):
    """
    t: time in days
    Mej: ejecta mass in solar masses
    
    Compute the radioactive heating rate L_in for a knova, as presented in:
    --> Korobkin+12, MNRAS 426, 1940-1949 (their equation 4)
    
    Uses constants T0_L_in = 1.3s, SIGMA_L_in = 0.11s defined in script.
    
    units: erg/s
    """
    
    l = (4e18)*Mej*Msun # assuming r-process mass = all ejecta
    l *= (0.5 - (1./np.pi)*np.arctan((t*86400-T0_L_in)/SIGMA_L_in))**1.3
    
    return l


def Barnes16_params(vej, Mej):
    """
    vej: ejecta velocity (as fraction of c)
    Mej: ejecta mass in solar masses
    
    Using the tabulated fit parameters of 
    --> Barnes+16, ApJ 829, 110 (their table 1)
    for the thermalization efficiency of a kilonova with a randomly-oriented
    magnetic field, linearly interpolate to obtain these parameters for the 
    input ejecta velocity and mass. 
    """
    ## interpolate Barnes+2016 to get thermalization fit parameters for KN    
    mej_barnes = np.log10(np.array([1e-3, 5e-3, 1e-2, 5e-2]))
    vej_barnes = np.array([0.1, 0.2, 0.3])    
    a_barnes = np.array([[2.01, 4.52, 8.16], [0.81, 1.90, 3.20], 
                         [0.56, 1.31, 2.19], [0.27, 0.55, 0.95]])
    b_barnes = np.array([[0.28, 0.62, 1.19], [0.19, 0.28, 0.45], 
                         [0.17, 0.21, 0.31], [0.10, 0.13, 0.15]])
    d_barnes = np.array([[1.12, 1.39, 1.52], [0.86, 1.21, 1.39], 
                         [0.74, 1.13, 1.32], [0.60, 0.90, 1.13]])
    
    fa_barnes = interp2d(vej_barnes, mej_barnes, a_barnes, kind='linear')
    fb_barnes = interp2d(vej_barnes, mej_barnes, b_barnes, kind='linear')
    fd_barnes = interp2d(vej_barnes, mej_barnes, d_barnes, kind='linear')
    
    a = fa_barnes(vej, np.log10(Mej))[0]
    b = fb_barnes(vej, np.log10(Mej))[0]
    d = fd_barnes(vej, np.log10(Mej))[0]
    
    return a, b, d
    

def eth(t, a, b, d):
    """
    t: time in days
    a, b, d: Barnes+2016 fit parameters (see Barnes16_params() above)
    
    computes the thermalization efficiency as a function of time, following
    --> Barnes+16, ApJ 829, 110 (their equation 36)
    
    units: unitless
    """
    
    return 0.36 * (np.exp(-a*t) + np.log(1. + 2.*b*(t**d))/(2.*b*(t**d)))


def td(kappa, mej, beta_ej):
    """
    kappa: (gray) characeristic ejecta opacity 
    mej: ejecta mass (solar masses)
    beta_ej: ejecta velocity v/c
    
    Compute diffusion time for a knova, as presented in 
    --> Chatzopoulos+12, ApJ 746, 121
    
    The constant 13.8 depends on the system geometry.
    
    units: s
    """
    return np.sqrt(2.*kappa*mej*Msun/(13.8*beta_ej*c*c))


def Lbol(t, eth, L_in, td):
    """
    t: time post-merger (any units as long as same as td)
    eth: thermalization efficiency
    L_in: input radioactive heating rate
    td: diffusion timescale 
    
    t, eth, and L_in **must** be arrays. td must be a scalar.
    
    Compute the bolometric luminosity for a kilonova under the assumptions of 
    (1) free (homologous) expansion, (2) negligible initial photosphere radius,
    and (3) negligible initial thermal energy following:
    --> Chatzopoulos+12, ApJ 746, 121 (their equation 3)
    
    units: erg/s
    """
    
    L_ret = []
    for i in range(1, len(t)+1):
        t_cut = t[:i]
        prefac = np.exp(-(t_cut[-1]/td)**2) # prefactor and integrand
        integrand = L_in[i-1]*eth[i-1]*np.exp((t_cut/td)**2) * (t_cut/td**2) 
        integ = trapz(integrand, x=t_cut)
        L_ret.append(prefac*integ)    
    return np.array(L_ret)


def Tphot(t, L_bol, beta_ej, Tc):
    """
    t: time in days
    L_bol: bolometric luminosity at this time
    vej: characteristic ejecta beta v/c
    Tc: critical temperature (i.e., temperature floor)
    
    Computes the temperature of the photosphere according to e.g. 
    --> Villar+2017, ApJ 851, L21 (their equation 4)
    under the assumption that the photosphere is well-described by a blackbody.
    """
    
    T = [max((L_bol[i]/(4.*np.pi*sigma_sb*(beta_ej*c*t[i]*86400)**2))**0.25, 
             Tc) for i in range(len(t))]
    return np.array(T)


def Rphot(t, L_bol, beta_ej, T, Tc):
    """
    t: time in days
    L_bol: bolometric luminosity at this time
    vej: characteristic ejecta beta v/c
    T: temperature at this time 
    Tc: critical temperature (i.e., temperature floor)
    
    Computes the radius of the photosphere according to e.g. 
    --> Villar+2017, ApJ 851, L21 (their equation 5)
    under the assumption that the photosphere is well-described by a blackbody.
    """
    R_ret = []
    for i in range(len(t)):
        if T[i] > Tc:
            R_ret.append(beta_ej*c*t[i]*86400)
        else:
            R_ret.append((L_bol[i]/(4.*np.pi*sigma_sb*Tc**4))**0.5)
            
    return np.array(R_ret)
    
    
 
               
               