#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 8 13:24:26 2020
@author: Nicholas Vieira
@lines.py

Load in line lists from the National Institute of Standards and Technology 
(NIST) Atomic Spectra Database (ASD), or, from the Kurucz & Bell line list. 

Currently, line lists from Kurucz need to be downloaded indpendently from the 
website https://www.cfa.harvard.edu/amp/ampdata/kurucz23/sekur.html 

An example of the format that should be used is shown on the github. 
"""

import numpy as np
import os
from astropy.table import Table
from astroquery.nist import Nist
import astropy.units as u

### NIST ASD ##################################################################

def query_lines_NIST(wavmin, wavmax, elements=None, fik_cut=0.05, write=None):
    """
    wavmin: lower wavelength bound in Angstroms
    wavmax: upper wavelength bound in Angstroms
    elements: string of elements of interest separated by a space character or 
              semi-colon, e.g. "Fe Ni Co"/"Fe;Ni;Co" OR a list of elements 
              e.g. ["Fe", "Ni", "Co"] 
              (optional; default None; which loads in "H I")
    fik_cut: minimum oscillator strength cut to keep a line (optional; default
             0.05)
    write: filename to write the table to (optional; default None, which will
           write no file)
    
    Submits an online query to obtain a table of lines from the NIST ASD for 
    the desired elements and wavelength ranges. Optionally, apply a minimum cut 
    on the oscillator strength to obtain only sufficiently strong lines. 
    Optionally writes the file after reading it in.
    """
    if type(elements) in [list, np.ndarray]:
        ecopy = elements.copy()
        elements = ""
        for e in ecopy: 
            elements += f"{e} "
    
    tab = Nist.query(wavmin*u.AA, wavmax*u.AA, elements, 
                     output_order="wavelength",
                     wavelength_type="vacuum")
    if len(tab) == 0:
        print("No lines found for given elements.")
        return 
    
    # mask out lines with no oscillator strength or no observed wavelength
    mask_nones = [type(t) != type(None) for t in tab["fik"].tolist()]
    tab = tab[mask_nones]
    mask_nones = [type(t) != type(None) for t in tab["Observed"].tolist()]
    tab = tab[mask_nones]
    
    # mask out lines with an alphabetic character in oscillator strength
    mask_numerical = []
    for t in tab:
        try:
            float(t["fik"])
            mask_numerical.append(True)
        except:
            mask_numerical.append(False)
    tab = tab[mask_numerical]
    
    # convert oscillator strength column to floats and apply cut to oscillator
    # strength
    tab["fik"] = [float(f) for f in tab["fik"].tolist()]
    tab = tab[tab["fik"] > fik_cut]

    # rename wavelength column 
    tab.rename_column("Observed","wavelength [AA]")
    
    if write:
        tab.write(write, format="ascii.csv", overwrite=True)
    
    return tab


def load_lines_NIST(table, wavmin=None, wavmax=None, fik_cut=None):
    """
    table: filename for table of NIST ASD lines to load in
    wavmin: lower wavelength bound in Angstroms (optional; default None)
    wavmax: upper wavelength bound in Angstroms (optional; default None) 
    fik_cut: minimum oscillator strength cut to keep a line (optional; default
             None)
    
    Load in a table of lines which was previously downloaded using the function
    query_lines_NIST() above. (In reality, can load in any line list as long 
    as the tables contain the columns 'wavelength [AA]' and 'fik', where the 
    latter is the oscillator strength).
    """
    
    tab = Table.read(table)
    
    # apply cut to wavelengths
    if wavmin:
        tab = tab[tab["wavelength [AA]"] > wavmin]
    if wavmax:
        tab = tab[tab["wavelength [AA]"] < wavmax]
    
    # apply cut to oscillator strengths
    if fik_cut:
        tab = tab[tab["fik"] > fik_cut]
        
    return tab

#### Kurucz line lists ########################################################
def fix_file_Kurucz(textfile, output=None):

    """
    textfile: textfile containing Kurucz line list (must be downloaded 
              separately by the user from 
              https://www.cfa.harvard.edu/amp/ampdata/kurucz23/sekur.html)\
    output: name for output fixed file (optional; default set below)
    
    Fixes some formatting inconsistencies which might be present in line lists
    downloaded from the above website so that astropy is not confused by the 
    formatting.
    """
        
    tf = open(textfile, "r")
    tf_rows = tf.readlines()
    tf_rows_list = [row.split() for row in tf_rows]
    tf.close()
    
    # check every row, make sure no empty columns
    # if an empty column is found, fill it with "-" 
    # need to do this so astropy doesn't get confused  
    for i in range(len(tf_rows_list)):
        # if following is true true, second col is log(gf)
        # --> missing wavelength in air column entry
        if float(tf_rows_list[i][1]) < 4.0: 
            tf_rows_list[i] = [tf_rows_list[i][0]]+["0.0"]+tf_rows_list[i][1:]
        # if still missing an entry
        if len(tf_rows_list[i]) == 11: 
            tf_rows_list[i] = tf_rows_list[i] + ["--"]          
  
        # add newline character to rows for writing 
        tf_rows_list[i] += ["\n"]
    
    if not(output):
        output = f"{os.path.splitext(textfile)[0]}_fixed.txt"
    
    newtf = open(output, "w+")
    newtf_rows = ["\t".join(t) for t in tf_rows_list]
    newtf.writelines(newtf_rows)
    newtf.close()
    print(f"new file written to {output}")
        
    
def load_lines_Kurucz(textfile, wavmin=None, wavmax=None):
    """
    textfile: textfile containing Kurucz line list (must be downloaded 
              separately by the user from 
              https://www.cfa.harvard.edu/amp/ampdata/kurucz23/sekur.html)
    wavmin: lower wavelength bound in Angstroms (optional; default None)
    wavmax: upper wavelength bound in Angstroms (optional; default None) 
    
    Load in a Kurucz line list from a file which was downloaded from the above 
    website. Automatically checks for formatting inconsistencies which can 
    confuse astropy, and fixes them, writing the fixed list to a new file.
    """
    # read in file
    try:
        tab = Table.read(textfile, format="ascii")
    except:
        print("empty columns --> fixing file")
        fix_file_Kurucz(textfile)
        tab = Table.read(f"{os.path.splitext(textfile)[0]}_fixed.txt", 
                         format="ascii")       
    
    # change wavelengths from nm to Angstroms 
    tab["col1"] *= 10
    
    # rename columns of table
    colnames = ["wavelength [AA]", 
                "wavelength_air [AA]", 
                "log_gf", "A [s^-1]", 
                "Elem.", "Element", "ionization", # e.g. 31.01, Ga, I
                "E_lower_lev [cm^-1]", "J lower", 
                "E_upper_lev [cm^-1]", "J upper", "Ref"]
    for i in range(len(tab.colnames)):
        tab.rename_column(f"col{i+1}", colnames[i])

    # apply cut to wavelengths
    if wavmin:
        tab = tab[tab["wavelength [AA]"] > wavmin]
    if wavmax:
        tab = tab[tab["wavelength [AA]"] < wavmax]

    return tab

#### find line following a given line in some line list #######################
def next_line(wavelen, line_table, offset=0):
    """
    wavelen: wavelength of interacting photon in angstroms
    line_table: Kurucz/NIST line list table
    offset: how many lines to look ahead (optional; default 0)
            (e.g., if you want not the next line but the line after, set 
            offset=1)
            
    Given either a NIST ASD line list, a Kurucz line list, or **any** line list
    which contains the column 'wavelength [AA]', and some wavelength for an 
    interacting photon, find the closest line in the list with a *longer* 
    wavelength than the photon. Optionally, provide the 2nd line, 3rd line, and
    so forth using the offset argument.
    """
    mask = line_table["wavelength [AA]"]>wavelen
    lines = line_table["wavelength [AA]"][mask].tolist()
    if len(lines) == 0:
        return 0
    argmin = np.argmin(np.array(lines) - wavelen)
    
    try:
        return lines[argmin+offset]
    except IndexError: # if no more lines...
        return 0

#### combine a line list from NIST ASD and Kurucz line list ###################
def combine_NIST_Kurucz(NIST_table_file, Kurucz_textfile, 
                        wavmin_NIST=None, wavmax_NIST=None,
                        wavmin_Kurucz=None, wavmax_Kurucz=None,
                        write=None):
    """
    NIST_table_file: filename for a NIST ASD line list
    Kurucz_textfile: filename for a Kurucz line list
    wavmin_NIST: minimum required wavelength for NIST ASD lines, in Angstroms
                 (optional; default None)
    wavmax_NIST: maximum allowed wavelength for NIST ASD lines, in Angstroms
                 (optional; default None)
    wavmin_Kurucz: minimum required wavelength for Kurucz lines, in Angstroms
                   (optional; defualt None)
    wavmax_Kurucz: maximum allowed wavelength for Kurucz lines, in Angstroms 
                   (optional; default None)
    write: whether to write the file or just load it in (optional; default 
           False, which does not write the file)
    
    Combine a NIST and Kurucz line list into a single file, for convenience. 
    """
    
    NIST_lines = load_lines_NIST(NIST_table_file, wavmin_NIST, wavmax_NIST)
    Kurucz_lines = load_lines_Kurucz(Kurucz_textfile, wavmin_Kurucz, 
                                     wavmax_Kurucz)
    
    # single-column table, just the wavelengths
    wavelens = NIST_lines["wavelength [AA]"].tolist()
    wavelens += Kurucz_lines["wavelength [AA]"].tolist()
    comb = Table(data=[wavelens], names=["wavelength [AA]"])
    
    if write:
        comb.write(write, format="ascii", overwrite=True)
    
    return comb

def load_lines_combined(line_file):
    """
    line_file: filename for combined NIST ASD + Kurucz line list
    
    Load in the combined line list.
    """
    return Table.read(line_file, format="ascii")
    
    

