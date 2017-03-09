# -*- coding: utf-8 -*-
"""
Configuration Settings
----------------------
Built for the VIRUS instrument as well as LRS2 on HET

@author: gregz
"""

config_dict = {'wvl_dict':'wl', 'specname':'sn', 'fsize':'fs', 'disp':'di',
               'ifucen_fn':'fn', 'cube_scale':'cs', 'fibmodel_bins':'bn',
               'wave_nbins':'wbn', 'default_fib':'dfn', 
               'fibmodel_sig': 'sig', 'fibmodel_pow':'pow', 
               'dark_mult':'dm', 'nfibers':'nf', 'cont_smooth':'contsmooth'}

# Bottom Amplifier for each side ---
Amps = ["LL","RU"]

# Connecting the bottom ampliefer with the top and total side
Amp_dict = {"LL": ["LU","L"], "RU": ["RL","R"]}

# Wavelength limits for each side as defined by the bottom amplifier id
virus_wl = {"LL": [3490,5500], "RU": [3490,5500]}
lrs2b_wl = {"LL": [3633,4655], "RU": [4550,7000]}
lrs2r_wl = {"LL": [6425,8440], "RU": [8230,10550]}
virusw_wl = {"LL": [4727,5503]}

# Dark multiplier for dark subtraction
virus_dm = {"LL": 1.0, "LU": 1.0, "RU": 1.0, "RL": 1.0}
lrs2b_dm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
lrs2r_dm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
virusw_dm ={"LL": 0.0, "LU": 0.0}

# Number of expected fibers
virus_nf = {"LL": 112, "LU": 112, "RU": 112, "RL": 112}
lrs2b_nf = {"LL": 140, "LU": 140, "RU": 140, "RL": 140}
lrs2r_nf = {"LL": 140, "LU": 140, "RU": 140, "RL": 140}
virusw_nf = {"LL": 133, "LU": 134}


# Name prefix for the normalized spectrum used for the wavelength solution
virus_sn = {"LL": "virus", "RU": "virus"}
lrs2b_sn = {"LL": "lrs2_uv", "RU": "lrs2_orange"}
lrs2r_sn = {"LL": "lrs2_red", "RU": "lrs2_farred"}
virusw_sn = {"LL": "virusw"}

# Name of the IFUcen file for fiber positions
virus_fn = {"LL": ["IFUcen_VIFU",30], "RU": ["IFUcen_VIFU",30]}
lrs2b_fn = {"LL": ["LRS2_B_UV_mapping.txt",4], "RU": ["LRS2_B_OR_mapping.txt",4]}
lrs2r_fn = {"LL": ["LRS2_R_NR_mapping.txt",4], "RU": ["LRS2_R_FR_mapping.txt",4]}
virusw_fn = {"LL": ["virusw_mapping.txt",4]}

# Pixel width in radius over which the fibermodel is defined
virus_fs = 8.
lrs2b_fs = 6.
lrs2r_fs = 6.
virusw_fs = 5.

# Number of bins for fibermodel fit
virus_bn = 15
lrs2b_bn = 11
lrs2r_bn = 11
virusw_bn = 7

# Number of bins for wavelength fit
virus_wbn = 21
lrs2b_wbn = 14
lrs2r_wbn = 14
virusw_wbn = 13

# Default fiber number for wavelength fit
virus_dfn = 42
lrs2b_dfn = 48
lrs2r_dfn = 48
virusw_dfn = 48

# The initial fibermodel sigma for initializing bins
virus_sig = 2.5
lrs2b_sig = 1.4
lrs2r_sig = 1.4
virusw_sig = 1.4

# The initial fibermodel power for initializing bins
virus_pow = 2.5
lrs2b_pow = 2.0
lrs2r_pow = 2.0
virusw_pow = 2.0

# Dispersion scale for making the Fe/CuFe files
virus_di = {"LL": 1.9, "RU": 1.9}
lrs2b_di = {"LL": 0.5, "RU": 1.2}
lrs2r_di = {"LL": 1.0, "RU": 1.0}
virusw_di = {"LL": 0.19}

# Cube pixel scale
virus_cs = 1.0
lrs2b_cs = 0.4
lrs2r_cs = 0.4
virusw_cs = 2.0

# Fraction of the image from which to select the default column for trace
frac = 0.47

# Continuum smoothing
virus_contsmooth = 25
lrs2b_contsmooth = 25
lrs2r_contsmooth = 25
virusw_contsmooth = 25