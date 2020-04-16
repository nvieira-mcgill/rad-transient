#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 02:06:29 2020
@author: Nicholas Vieira
@collect_results.py

Produce a kilonova light curve using the results from previously completed runs
of the Monte Carlo code rad_transient.py.
"""

import numpy as np
import glob
import sys

import matplotlib.pyplot as plt 
from matplotlib import rc
rc('text', usetex=True)

## collect the data ###########################################################
try:
    DATADIR = str(sys.argv[1])
except IndexError:
    DATADIR = "./previous_runs"
allfiles = glob.glob(f"{DATADIR}/*results*.npy")
t_arrival_all = []
final_e_all = []
final_freq_all = []

for a in allfiles:
    arr = np.load(a)
    t_arrival_all += arr[0].tolist()
    final_e_all += arr[1].tolist()
    final_freq_all += arr[2].tolist()
    del arr

## make a "light curve" from the Monte Carlo ##################################
tmin = np.floor(min(t_arrival_all)) # start
tmax = 15

# mask results which escape after 15 days
mask = np.array(t_arrival_all) < tmax
t_arrival_all = np.array(t_arrival_all)[mask]
final_e_all = np.array(final_e_all)[mask]
final_freq_all = np.array(final_freq_all)[mask]
dt = 0.5 # 0.5 days
t_lc = np.linspace(tmin, tmax, int((tmax-tmin)/dt)+1)
t_lc = np.logspace(np.log10(tmin), np.log10(tmax), int((tmax-tmin)/dt)+1)
lightcurve = np.zeros(t_lc.shape)

t_arrival_working = t_arrival_all.copy()
for i in range(len(t_lc)):
    t = t_lc[i]
    for j in range(len(t_arrival_all)):
        if t <= t_arrival_all[j] < t+dt:
            lightcurve[i] += final_e_all[j]/(dt*86400)

# divide by number of packets 
lightcurve /= len(final_e_all)

# subtract light crossing time from times array
t_lc -= 1e16/(3e10) / 86400

## make a prediction using the model for Lbol for comparison ##################
import heating

## ejecta parameters
M_EJ = 0.02 # total ejecta mass [solar masses]
BETA_EJ = 0.2 # ejecta characteristic velocity v/c
KAPPA_EJ = 1.0 # ejecta characteristic wavelength-averaged opacity [cm^2 g^-1]
YE = 0.4 # electron fraction (exact value arbitrary, only matters if above or
         # below 0.25)
TFLOOR = 2500 # temperature floor [K]

t_arr = np.logspace(np.log10(min(t_lc)), np.log10(max(t_lc)), 100)

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

## plotting the light curves ##################################################
fig = plt.figure(figsize=(26,16))

# the MC light curve
Nph_str = f"{len(t_arrival_all):.1e}" # nice string formatting
Nph_str = Nph_str[0] + r"{\times}10^{" + Nph_str[-1] +"}"
mask = lightcurve > 0
plt.step((t_lc + 0.5*dt)[mask], lightcurve[mask], marker="", lw=6, ls="-", 
         color="#0485d1", zorder=2, 
         label=r"Monte Carlo [$N_{\mathrm{ph}} = "+Nph_str+"$]")
         
# the model light curve 
plt.plot(t_arr, Lbol_arr, marker="", lw=6, ls="--", alpha=0.8, 
         color="#cb416b", zorder=1, label="Equation (3)")

plt.xlabel(r"$t$ [days]", fontsize=30)
plt.ylabel(r"Luminosity [erg $\mathrm{s^{-1}}$]", fontsize=30)
plt.yscale("log")
plt.xticks(size=35)
plt.yticks(size=35)
plt.xlim(min(t_arr), max(t_lc))
plt.ylim(min(lightcurve[lightcurve>0]), 5e40)
plt.grid(True, which='major', linewidth=1.8)
plt.grid(True, which='minor', linewidth=0.7)
plt.legend(loc='best', fontsize=35)
plt.savefig("results_lightcurve.pdf", bbox_inches="tight")