#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 01:42:47 2020
@author: Nicholas Vieira
@collect_MC_info.py

Collect Monte Carlo diagnostics such as computation times from previously 
completed runs produced by rad_transient.py.
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
allfiles = sorted(glob.glob(f"{DATADIR}/*MC_info*.npy"))
ppack_files = sorted(glob.glob(f"{DATADIR}/*MC_info_ppack*.npy"))
singlepack_files = [a for a in allfiles if not(a in ppack_files)]

### single-packet files #######################################################
## should all have the same NSHELL, NTSTEP, and NPACK for valid comparison
# assumed that they do
MC_time_all = []
arr = np.load(singlepack_files[0]) # load in first file
NSHELL = arr[1] # take NSHELL, NTSTEP and NPACK from first file
NTSTEP = arr[2] # assuming they are the same for all files
NPACK = arr[3]
del arr

for s in singlepack_files:
    arr = np.load(s)
    MC_time_all.append(arr[0])
    del arr

## plotting histograms 
# total MC time
fig = plt.figure(figsize=(16,14)) 
plt.hist(MC_time_all, 100, color="#0165fc", histtype="step", lw=4)
plt.xlabel("Monte Carlo time for all packets [s]", fontsize=32)
plt.ylabel("Counts", fontsize=32)
plt.xticks(size=35)
plt.yticks(size=35)
plt.grid()
plt.savefig(f"MC_time_all.pdf", bbox_inches="tight")
   
### per-packet files ##########################################################
MC_time_ppack_success_all = [] 
MC_time_ppack_failed_all = []

for p in ppack_files:
    arr = np.load(p)
    try: # sometimes an array, sometimes a list - not sure why 
        MC_time_ppack_success_all += arr[0].tolist()
    except AttributeError: # if a list, will have no method .tolist()
        MC_time_ppack_success_all += arr[0]
    
    try:
        MC_time_ppack_failed_all += arr[1].tolist()
    except AttributeError: # if a list, will have no method .tolist()
        MC_time_ppack_failed_all += arr[1]
    del arr

# compute success rate
success_rate = len(MC_time_ppack_success_all)/(len(MC_time_ppack_success_all)+
                   len(MC_time_ppack_failed_all))

## plotting histograms 
fig = plt.figure(figsize=(16,14)) 
# time taken for successful escapes 
hbins = np.logspace(np.log10(np.nanmin(MC_time_ppack_success_all)), 
                    np.log10(np.nanmax(MC_time_ppack_success_all)), 
                    100)
plt.hist(MC_time_ppack_success_all, bins=hbins, 
         color="#0504aa", histtype="step", lw=4, 
         label=f"Escaped [{success_rate*100:.2f}"+"\%]")
# for failed escapes
hbins = np.logspace(np.log10(np.nanmin(MC_time_ppack_failed_all)), 
                    np.log10(np.nanmax(MC_time_ppack_failed_all)), 
                    100)
plt.hist(MC_time_ppack_failed_all, bins=hbins, 
         color="#ff028d", histtype="step", lw=4,
         label=f"Lost/absorbed [{(1-success_rate)*100:.2f}"+"\%]")
         
plt.xlabel("Monte Carlo time for one packet [s]", fontsize=32)
plt.ylabel("Counts", fontsize=32)
plt.xscale("log")
plt.xticks(size=35)
plt.yticks(size=35)
plt.grid(True, which='major', linewidth=1.8)
plt.grid(True, which='minor', linewidth=0.7)
plt.legend(loc='best', fontsize=35)
plt.savefig("MC_time_ppack_all.png", bbox_inches="tight")

### compare arrival times with computation times for individual packets #######
results_files = sorted(glob.glob(f"{DATADIR}/*results*.npy"))
t_arrival_all = []
for r in results_files:
    arr = np.load(r)
    t_arrival_all += arr[0].tolist()
    del arr
    
# get the arrival time percentiles
qs = [50,90]
pers = np.percentile(t_arrival_all, qs)
lines = ["--", "-"]

# plot
fig = plt.figure(figsize=(16,14)) 
plt.plot(t_arrival_all, MC_time_ppack_success_all, ls="", marker="o",
         ms=2, color="#0504aa")
for i in range(len(qs)):
    plt.axvline(pers[i], lw=3, color="k", label=f"{qs[i]}\%", ls=lines[i])
    
plt.xlabel("physical arrival time [days]", fontsize=32)
plt.ylabel("Monte Carlo time for one packet [s]", fontsize=32)
plt.xticks(size=35)
plt.yticks(size=35)
plt.xlim(left=min(t_arrival_all), right=15)
plt.ylim(bottom=-0.1)
plt.grid(True, which='major', linewidth=1.8)
plt.grid(True, which='minor', linewidth=0.7)
plt.legend(loc='best', fontsize=35)
plt.savefig("MC_time_ppack_tarrival_compare.pdf", bbox_inches="tight")