# -*- coding: utf-8 -*-
"""
Configuration Settings
----------------------
Built for the VIRUS instrument as well as LRS2 on HET

@author: gregz
"""

# Bottom Amplifier for each side ---
Amps = ["LL","RU"]

# Connecting the bottom ampliefer with the top and total side
Amp_dict = {"LL": ["LU","L"], "RU": ["RL","R"]}

# Wavelength limits for each side as defined by the bottom amplifier id
virus_wl = {"LL": [3490,5500], "RU": [3490,5500]}
lrs2b_wl = {"LL": [3633,4655], "RU": [4550,7000]}
lrs2r_wl = {"LL": [3490,5500], "RU": [3490,5500]}

# Dark multiplier for dark subtraction
virus_dm = {"LL": 1.0, "LU": 1.0, "RU": 1.0, "RL": 1.0}
lrs2b_dm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
lrs2r_dm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}

# Number of expected fibers
virus_nf = {"LL": 112, "LU": 112, "RU": 112, "RL": 112}
lrs2b_nf = {"LL": 140, "LU": 140, "RU": 140, "RL": 140}
lrs2r_nf = {"LL": 140, "LU": 140, "RU": 140, "RL": 140}

# Name prefix for the normalized spectrum used for the wavelength solution
virus_sn = {"LL": "virus", "RU": "virus"}
lrs2b_sn = {"LL": "lrs2_uv", "RU": "lrs2_orange"}
lrs2r_sn = {"LL": "lrs2_red", "RU": "lrs2_farred"}

# Name of the IFUcen file for fiber positions
virus_fn = {"LL": ["IFUcen_VIFU",30], "RU": ["IFUcen_VIFU",30]}
lrs2b_fn = {"LL": ["LRS2_B_UV_mapping.txt",4], "RU": ["LRS2_B_OR_mapping.txt",4]}
lrs2r_fn = {"LL": "LRS2_R_NR_mapping.txt", "RU": "LRS2_R_FR_mapping.txt"}

# Pixel width in radius over which the fibermodel is defined
virus_fs = 8.
lrs2b_fs = 6.
lrs2r_fs = 6.

# Number of bins for fibermodel fit
virus_bn = 15
lrs2b_bn = 11
lrs2r_bn = 11

# The initial fibermodel sigma for initializing bins
virus_sig = 2.5
lrs2b_sig = 1.4
lrs2r_sig = 1.4

# The initial fibermodel power for initializing bins
virus_pow = 2.5
lrs2b_pow = 2.0
lrs2r_pow = 2.0

# Dispersion scale for making the Fe/CuFe files
virus_di = {"LL": 1.9, "RU": 1.9}
lrs2b_di = {"LL": 0.5, "RU": 1.2}
lrs2r_di = {"LL": 1.2, "RU": 1.2}

# Cube pixel scale
virus_cs = 1.0
lrs2b_cs = 0.4
lrs2r_cs = 0.4

# Fraction of the image from which to select the default column for trace
frac = 0.47

# Continuum smoothing
virus_contsmooth = 25
lrs2b_contsmooth = 25
lrs2r_contsmooth = 25
