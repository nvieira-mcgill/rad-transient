#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 18:52:31 2020
@author: Nicholas Vieira
@rad_transient.py

Monte Carlo code to produce a light curve/spectra for a radioactivel-powered
transient such as a supernova or kilonova. Based primarily on the work of:

--> Bulla 2019, MNRAS 489, 5037-5045
--> Bulla et al. 2015, MNRAS 450, 967–981
--> Kasen et al. 2006, ApJ 651, 366-380
--> Mazzali & Lucy 1993, A&A 279, 447-456
--> Tanaka et al. 2019, arXiv:1906.08914

See the accompanying report on the github 
(https://github.com/nvieira-mcgill/rad-transient/) for details.

"""

## imports
# basic
from timeit import default_timer as timer
import time
import numpy as np

# import relevant constants from astropy, convert to cgs
from astropy.constants import h, c, M_sun
h = h.cgs.value # Planck's constant
c = c.cgs.value # speed of light
M_sun = M_sun.cgs.value # mass of the sun

# my modules
import lines # querying, loading in line lists
import profiles # density profiles 
import distribs # for thermal emissivity distribution
import heating # computing heating rates/thermalization/luminosity/photospheres
import event_selection # selecting an event during the Monte Carlo
import electron_scattering # scatter a photon off of an electron
import line_scattering # scatter a photon in resonance with an atomic line

#### PARAMETERS ###############################################################
## ejecta parameters
M_EJ = 0.02 # total ejecta mass [solar masses]
BETA_EJ = 0.2 # ejecta characteristic velocity v/c
KAPPA_EJ = 1.0 # ejecta characteristic wavelength-averaged opacity [cm^2 g^-1]
YE = 0.4 # electron fraction (exact value arbitrary, only matters if above or
         # below 0.25)
TFLOOR = 2500 # temperature floor [K]
A = 1e35 # scaling for density profile 

## Monte Carlo parameters
T_INIT = 0.1 # start time [days]
TMAX = 15 # maximum allowed time for photon to escape [days]
RMAX = 1e16 # edge of computational grid, i.e. outermost shell [cm]
NTSTEP = 10 # number of time steps
NSHELL = 100 # number of radial shells
NPACK = 1 # no. of photon packets to initialize at each time step
TIMESLEEP = 0#5 # time [s] to sleep between packet escapes

#### LINE LISTS ###############################################################
WAVMIN = 0 # min wavelength to acquire from line lists [Angstroms]
WAVMAX = 25000 # max wavelength to acquire from line lists [Angstroms]

## elements to query from NIST
#ELEMENTS = ["Ge", #32 
#            "Rb", #37
#            "Sr", #38
#            "Y", #39
#            "Zr", #40
#            "Nb", #41
#            "Mo", #42 
#            "Tc", #43
#            "Ru", #44
#            "Rh", #45
#            "Pd", #46
#            "Ag", #47
#            "Cs", #55
#            "Ba", #56
#            "Ir", #76
#            "Os", #77
#            "Pt"] #78


## get the NIST line list and combine it with the Kurucz line list,
## using NIST for < 2200 AA and Kurucz for 2200-25000 AA
#line_table_NIST = lines.query_lines_NIST(WAVMIN, WAVMAX, ELEMENTS,
#                                         write="linelists/NIST_rprocess.csv")
#line_table = lines.combine_NIST_Kurucz("linelists/NIST_rprocess.csv",
#                                "linelists/rprocess_nohdr_bigrange_fixed.txt",
#                                wavmin_NIST=WAVMIN,
#                                wavmax_NIST=2200,
#                                wavmin_Kurucz=2200,
#                                wavmax_Kurucz=WAVMAX,
#                                write="NIST_and_Kurucz_rproc.csv")

# or, just load in the combined table
line_table = lines.load_lines_combined("NIST_and_Kurucz_rproc.csv")

#### RADIUS AND TEMPERATURE OF PHOTOSPHERE AT EACH TIME STEP ##################
# create <NTSTEP> logarithmically-spaced time steps
t_arr = np.logspace(np.log10(T_INIT), np.log10(TMAX), NTSTEP+1)

# compute radioactive heating rate over time
Lin_arr = heating.Lin(t_arr, M_EJ)

# get thermalization efficiency fit parameters, compute efficiency over time
a, b, d = heating.Barnes16_params(BETA_EJ, M_EJ)
etherm_arr = heating.eth(t_arr, a, b, d)

# compute diffusion timescale
td = heating.td(KAPPA_EJ, M_EJ, BETA_EJ)

# compute bolometric luminosity over time 
Lbol_arr = heating.Lbol(t_arr, etherm_arr, Lin_arr, td/86400.)

# discard first value in all of the arrays to avoid singularity (Lbol(0) = 0) 
t_arr = t_arr[1:]
Lin_arr = Lin_arr[1:]
etherm_arr = etherm_arr[1:]
Lbol_arr = Lbol_arr[1:]

# compute photosphere temperature and radius over time 
T_phot_arr = heating.Tphot(t_arr, Lbol_arr, BETA_EJ, TFLOOR)
R_phot_arr = heating.Rphot(t_arr, Lbol_arr, BETA_EJ, T_phot_arr, TFLOOR)

# get **initial** density as a function of r 
rhoi_arr = profiles.rho_init(R_phot_arr, a=A)

# create <NSHELL> logarithmically-spaced concentric shells
r_ejecta = np.logspace(np.log10(R_phot_arr[0]), np.log10(RMAX), NSHELL)

#### MISC. ####################################################################
unif_z = np.linspace(0, 0.99999) # uniform dist from [0, 1)

#### MONTE CARLO ##############################################################
"""
At each time step, using the bolometric luminosity, temperature, and radius of 
the photosphere, <NPACK> photon packets are generated at the surface of the 
photosphere. Their initial propagation direction is isotropically away from
the photosphere.

These packets are followed in both time and space until they reach the 
computational boundary, at which point they are collected. Their energy, the 
frequency of the (identical) photons in the packet, and the arrival time are 
recorded.
"""

## record the time taken for the entire Monte Carlo
start = timer()

## arrays to be filled
# results
t_arrival = [] # time when photon exits computational boundary
final_e = [] # energy upon exit
final_freq = [] # frequency upon exit
final_mu = [] # cosine of propagation angle upon exit 

# MC info
MC_time_ppack_success = [] # time taken for each packet which escapes
MC_time_ppack_failed = [] # time taken for each packet which does not escape

## time steps #################################################################
for j in range(1, len(t_arr)-1): # at each time step  
    ## basics    
    t = t_arr[j] # time  
    T_phot = T_phot_arr[j] # temperature of photosphere
    R_phot = R_phot_arr[j] # radius of photosphere
    
    # compute the energy injected in the last time step
    # don't forget to convert days to seconds 
    E_inject = Lin_arr[j]*etherm_arr[j]*(t_arr[j]-t_arr[j-1])*86400 
    print(f"{E_inject:.1e} erg injected")
    
    E_pack = E_inject/NPACK # compute energy per photon packet
    print(f"{E_pack:.1e} erg per packet\n")
    
    ## setup the photon packet(s)
    # randomly sample to obtain the frequency of the photons in each packet
    freq_packs = distribs.gen_thermal_emissivity(T_phot, t, YE, n=NPACK)/h

    # compute no. of photons in each packet 
    nph_packs = E_inject/(h*freq_packs)
    
    # compute the initial prop. direction of the packet(s) (isotropic outwards)
    mu_packs = np.sqrt(1 - np.random.choice(unif_z, size=NPACK))
        
    ## for each packet ########################################################
    for n in range(NPACK):
        start_pack = timer() # also record the time taken for each packet alone        
        # setup the packet
        r_n = R_phot # starting position (photosphere surface)
        t_n = t # starting time
        e_n = E_pack # energy of packet (same for all packets)
        freq_n = freq_packs[n] # frequency       
        wavelen_n = c/freq_n * 1e4 # wavelength (in microns)
        mu_n = mu_packs[n] # cosine of propgation angle
        nph_n = E_pack/(h*freq_packs) # number of photons in packet
        
        ### propagate packet in both space and time until it reaches <RMAX> ###
        while r_n < RMAX:
            ## find current shell (coordinate i in r-grid) 
            ## find current time step (coordinate j in time-grid)
            idx_current = len(r_ejecta) - len(r_ejecta[r_ejecta >= r_n])
            idx_time = len(t_arr) - len(t_arr[t_arr >= t_n])
            
            if idx_current == len(r_ejecta):
                print("** packet was back-scattered beyond initial "+
                      "photosphere --> packet absorbed and lost\n")
                end_pack = timer() # record time taken
                MC_time_ppack_failed.append(end_pack-start_pack)
                break
            if idx_time == len(t_arr):
                print("** escape time limit reached by packet --> packet "+
                      "absorbed and lost\n")
                end_pack = timer() # record time taken
                MC_time_ppack_failed.append(end_pack-start_pack)
                break
            
            print(f"i = {idx_current}/{len(r_ejecta)}\t\t"+
                  f"j = {idx_time}/{len(t_arr)}")

            # compute distance to current shell + nearest shell **ahead** 
            r_current = r_ejecta[r_ejecta <= r_n][-1] # current shell (i)
            r_ahead = r_ejecta[r_ejecta > r_n][0] # shell ahead (i+1)
            drad_to_sh_ahead = r_ahead - r_n
            print(f"t = {t_n:.2e} days\tr = {r_n:.2e} cm")
            print(f"r_i = {r_current:.2e} cm\t"+
                  f"r_i+1 = {r_ahead:.2e} cm\t"+
                  f"d_r_i+1 = {drad_to_sh_ahead:.2e} cm")
            
            # get density and temperature at current r and time  
            rho_n = profiles.rho(r_n, t_n, t0=T_INIT, rho_i=rhoi_arr)[idx_time]
            T_n = T_phot_arr[idx_time]
            print(f"rho_ij = {rho_n:.2e} g cm**-3\t\t"+
                  f"T_j = {T_n:.1f} K")
            
            ## select the event and compute distance travelled  ###############
            mu_old = mu_n # record incident angle 
            v_shell = BETA_EJ*c # needs to be changed????
            if mu_n < 0: # if back-scattering
                dist_to_sh_next = (r_n - r_current)/abs(mu_n)
            else:
                dist_to_sh_next = drad_to_sh_ahead/mu_n
            event, dtravel = event_selection.event_select(t_n, wavelen_n, 
                                                          rho_n, v_shell/c,
                                                          mu_n, 
                                                          dist_to_sh_next,
                                                          line_table, 
                                                          YE)
            if event == "electron scatter": 
                print("electron scattering")
                # update direction, energy
                mu_n, e_n  = electron_scattering.scatter(v_shell/c, mu_n, e_n)
                print(f"mu_new = {mu_n:.2f}")
                
            ## if line scattering
            elif type(event) == float:
                print(f"interacting with {event:.2f}A line")
                # compute optical depth of the line in Sobolev approximation
                tau_line = event_selection.tau_line_sob(t_n, wavelen_n, rho_n, 
                                                        event/1e4, 
                                                        dtravel, YE)
                # update direction, energy
                mu_n, e_n = line_scattering.scatter(tau_line, v_shell/c, 
                                                    mu_n, e_n)
                if mu_n == mu_old: # if did not interact (low probability)
                    print("no scatter!")
                    print("exiting shell")
                    dtravel = dist_to_sh_next # go to next shell instead
                else:
                    print(f"mu_new = {mu_n:.2f}")    
                    
            ## if exiting the shell
            else:     
                print("exiting shell")
            
            # travel a radial distance <dtravel * mu_n>
            # can be forward or backward
            print(f"distance travelled = {dtravel:.1e} cm")
            print(f"radial distance travelled = {dtravel*mu_n:.1e} cm\n\n")
            r_n += dtravel*mu_n # travel a radial distance <dtravel * mu_n>
            t_n += (abs(dtravel)/c) / 86400 # elapse some time dtravel/c
                
        ## when the packet has escaped
        if r_n >= RMAX:
            t_arrival.append(t_n)
            final_e.append(e_n)
            final_freq.append(e_n/(nph_n*h))
            final_mu.append(mu_n)
            print(f"** packet escaped at t = {t_n:.1e} days!\n")
            end_pack = timer()
            MC_time_ppack_success.append(end_pack-start_pack)
        time.sleep(TIMESLEEP) # pause before moving on to next photon packet       


#### RECORD RESULTS ###########################################################
# actual physical results
final_freq = [f[0] for f in final_freq]
results = np.array([t_arrival, 
                    final_e, 
                    final_freq, 
                    final_mu])
# Monte Carlo diagnostics
end = timer()
MC_info = np.array([end-start, # total MC time
                    NSHELL,
                    NTSTEP,
                    NPACK])
MC_info_ppack = np.array([MC_time_ppack_success, # times per packet
                          MC_time_ppack_failed])

print(f"\nTIME FOR ENTIRE MONTE CARLO = {end-start:.3f} s")

## save files
from astropy.time import Time
now = Time.now().isot.replace(":","-")
np.save(f"rad_transient_results_{now}.npy", results) 
np.save(f"rad_transient_MC_info_{now}.npy", MC_info)
np.save(f"rad_transient_MC_info_ppack_{now}.npy", MC_info_ppack)

    
