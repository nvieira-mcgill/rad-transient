# rad-transient
Basic Monte Carlo code for simulating a radioactively-powered astrophysical transient such as a supernova or kilonova. The only radiative processes considered are electron scattering and bound-bound interactions. Written for McGill course PHYS642 - Radiative Processes in Astrophysics. The methods used are based primarily on the work of:

* Bulla 2019, MNRAS 489, 5037-5045 https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.5037B/abstract
* Bulla et al. 2015, MNRAS 450, 967–981
* Kasen et al. 2006, ApJ 651, 366-380
* Mazzali & Lucy 1993, A&A 279, 447-456
* Tanaka et al. 2019, arXiv:1906.08914

For more details, see the accompanying report in the file `PHYS642_Project2.pdf`


**main script**: `rad_transient.py`

This script loads in a line list containing an assortment of <em>r</em>-process elements and runs a Monte Carlo to propagate packets of photons through an ejecta typical of a blue, lanthanide-poor kilonova. 

**Other important scripts:**

