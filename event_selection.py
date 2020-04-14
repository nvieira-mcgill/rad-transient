#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:40:25 2020
@author: Nicholas Vieira
@event_selection.py

Given some photon in some supernova/kilonova ejecta, select the interaction
which the photon undergoes, where options are:
    (1) electron scattering 
    (2) scattering off of a line 
    (3) exiting the shell (distance interval) of the grid
    
The events are chosen using the schematic originally outlined in:
--> Mazzali & Lucy 93, A&A 279, 447-456 (their section 3.4)
"""

import numpy as np
import opac
import lines
import doppler

from astropy.constants import c
c = c.cgs.value

#### DISTRIBUTIONS ############################################################
z_dist = np.linspace(0,0.999,100)

#### OPTICAL DEPTHS ###########################################################
def tau_event():
    """
    Generate a random 'event' optical depth = -ln(z), where z lies in the 
    range [0, 1)
    """
    return -np.log(1 - np.random.choice(z_dist))


def tau_elec(t, rho, dist):
    """
    t: time post-merger in days
    rho: density in g cm**-3
    dist: distance to the interaction in cm
    
    Compute the optical depth to electron scattering for some known distance
    to the interaction.
    """
    return opac.kappa_es(t) * rho * dist


def tau_line_sob(t, wavelen, rho, line, dist, Ye):
    """
    t: time post-merger in days
    wavelen: rest-frame wavelength of the photon in microns
    rho: density in g cm**-3
    line: rest wavelength of the line of interest, in microns
    dist: distance to the interaction in cm
    
    Compute the optical depth of a specific line under the Sobolev 
    approximation.
    """
    if Ye < 0.25: # red 
        kappa_bb = opac.kappa_bb_lanthrich(t, 
                                      prefac=opac.prefac_bb_lanthrich(wavelen))
    else: # blue
        kappa_bb = opac.kappa_bb_lanthpoor(t, 
                                      prefac=opac.prefac_bb_lanthpoor(wavelen))
    return kappa_bb * rho * dist


#### EVENT SELECTION SCHEME ###################################################
def event_select(t, wavelen, rho, beta, mu, dist_to_shell, line_table, Ye):
    """
    t: time post-merger in days
    wavelen: rest-frame wavelength of photon in **microns**
    rho: density in g cm**-3
    beta: ejecta velocity v/c in shell
    mu: cosine of incident angle of photon 
    dist_to_shell: distance which must be travelled to exit the shell
    line_table: table of lines of interest (must have, at minimum, a column 
                'wavelen [AA]' with the wavelengths in angstroms)
    Ye: electron fraction
    
    Select an event to occur (electron scattering, line scattering, or 
    no interaction) using the scheme originally outlined in 
    --> Mazzali & Lucy, A&A 279, 447-456 (their section 3.4)
    """
    
    tau_r = tau_event() # draw an "event" optical depth
    
    # compute distance to electron scattering 
    dist_e = tau_r/(opac.kappa_es(t)*rho)   
    
    ## find next line interaction, compute distance
    # first, find wavelength of photon in co-moving frame (in microns)
    wavelen_cmf = doppler.wavelen_cmf(wavelen, beta, mu)
    print(f"{wavelen*1e4:.2f} AA\t{wavelen_cmf*1e4:.2f} AA"+
          f"\t{(wavelen_cmf - wavelen)*1e4:.2f} AA")
    # find next line in the line list (don't forget to convert to angstroms)
    nl = lines.next_line(wavelen_cmf*1e4, line_table) 
    print(f"next line: {nl:.2f} AA")
    # finally, compute distance to next line
    v_l = c*(nl - wavelen_cmf*1e4)/nl # distance in v-space
    dist_l = v_l*t*86400 # actual distance
    
    distances = [dist_e, dist_l, dist_to_shell]
    print(f"{distances[0]:.2e} cm\t{distances[1]:.2e} cm"+
          f"\t{distances[2]:.2e} cm")
    
    tau_e = tau_elec(t, rho, dist_l)
    tau_l = tau_line_sob(t, wavelen, rho, nl, dist_l, Ye)
    print(f"{tau_e:.2e}\t{tau_l:.2e}\t{tau_r:.2e}\n")
    
    ## finally, select the event  
    # if distance to electron scattering is the shortest
    if np.argmin(distances) == 0: 
        return "electron scatter", dist_e # electron scattering
    
    # else if distance to line scattering is the shortest
    elif np.argmin(distances) == 1: 
        # if sum of scattering optical depths does not exceed the event optical 
        # depth, go on to next line
        skip=1 # how many lines to skip ahead 
        while tau_e + tau_l < tau_r:
            # does not exceed event optical depth --> go to next line
            nl = lines.next_line(wavelen_cmf*1e4, line_table, offset=skip)
            
            if nl == 0: # returned when there are no valid lines
                return "exit shell", dist_to_shell
            
            # compute distance to this next line
            v_l = c*(nl - wavelen_cmf*1e4)/nl
            dist_l = v_l*t*86400
            
            # if the distance to this new line exceeds the distance to the 
            # shell, exit the shell
            if dist_l > dist_to_shell: 
                return "exit shell", dist_to_shell # exit the shell
            
            # otherwise, recompute optical depths, and try again
            tau_e = tau_elec(t, rho, dist_l)
            tau_l = tau_line_sob(t, wavelen, rho, nl, dist_l, Ye)
            skip += 1 # move another line ahead
        return nl, dist_l
        
    else: # if dist_to_shell is shortest 
        return "exit shell", dist_to_shell
