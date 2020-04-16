# rad-transient
Basic Monte Carlo code for simulating a radioactively-powered astrophysical transient such as a supernova or kilonova. The only radiative processes considered are electron scattering and bound-bound interactions. Written for the McGill course PHYS642 - Radiative Processes in Astrophysics. The methods used are based primarily on the work of:

* [Bulla 2019, MNRAS 489, 5037-5045](https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.5037B/abstract)
* [Bulla et al. 2015, MNRAS 450, 967–981](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..967B/abstract)
* [Kasen et al. 2006, ApJ 651, 366-380](https://ui.adsabs.harvard.edu/abs/2006ApJ...651..366K/abstract)
* [Mazzali & Lucy 1993, A&A 279, 447-456](https://ui.adsabs.harvard.edu/abs/1993A%26A...279..447M/abstract)
* [Tanaka et al. 2019, arXiv:1906.08914](https://ui.adsabs.harvard.edu/abs/2019arXiv190608914T/abstract)

For more details, see the accompanying report in the file `PHYS642_Project2.pdf`

**main script**: `rad_transient.py`

This script loads in a line list containing an assortment of <em>r</em>-process elements and runs a Monte Carlo to propagate packets of photons through an ejecta typical of a blue, lanthanide-poor kilonova. This line list is already included in `NIST_and_Kurucz_rproc.csv`. Lines are acquired from a combination of the National Institute of Standards and Technology (NIST) [Atomic Spectra Database (ASD)](https://physics.nist.gov/PhysRefData/ASD/lines_form.html) and the [Kurucz & Bell](https://www.cfa.harvard.edu/amp/ampdata/kurucz23/sekur.html) line list.

**Other important scripts:**

`collect_results.py` and `collect_MC_info.py` can be run from the command line, where the argument supplied should be a directory containing files output by `rad_transient.py`, e.g. 

`>>> python3 collect_results.py ~/my_data_directory/` 

If no argument is supplied, the script will look for files in a directory `./previous_runs`. This directory is included in a zip file on this github, and includes previous runs which amount to around 30,000 photon packets. 

Running the first of these scripts (`collect_results.py`) will produce a light curve from the runs (Figure 5 of the accompanying report) and running the second (`collect_MC_info.py`) will plot a series of diagnostics on e.g. computation times for the Monte Carlo itself (Figures 6-8). 
