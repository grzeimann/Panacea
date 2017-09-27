# -*- coding: utf-8 -*-
"""
Configuration Settings
----------------------
Built for the VIRUS instrument as well as LRS2 on HET
Works for VIRUS-W as well.

@author: gregz
"""

import numpy as np

# Output Directory
output = "reductions"
limited_output_files = True

# Common Options
adjust_trace = True # Adjust science trace for shifts
trace_from_sci = False # if adjust_trace is true, this is overridden
use_trace_ref = True # Use default fiber files to always recover desired fibers
refit_fiber_to_fiber = False # Adjust fiber to fiber using science
adjust_ftf = False # Adjust fiber to fiber using look-up table
use_other_sky = False # Use another sky background for desired sky subtraction
use_pixelflat = True
check_fibermodel = True
check_wave = True
cosmic_iterations = 1
make_model_image = False

# OPERATIONS
twi_prepare_image = True
twi_get_trace = True
twi_subtract_background = False
twi_get_fibermodel = True
twi_get_wavelength_solution = True
twi_get_fiber_to_fiber = True

sci_subtract_background = False
sci_sky_subtraction = True
sci_remeasure_fibermodel = False
sci_significance_map = True

# Configuration Directories
#rootdir = "/work/03946/hetdex/maverick"
#virusconfig = "/work/03946/hetdex/maverick/virus_config"
#darkpath = "/work/03730/gregz/maverick/lib_dark/march"
#biaspath = "/work/03730/gregz/maverick/lib_bias/march"
#pixflatpath = "/work/03946/hetdex/maverick/virus_config/PixelFlats"


rootdir = "/Users/gregz/cure/virus_raw"
virusconfig = "/Users/gregz/cure/virus_early/virus_config"
darkpath = "/Users/gregz/cure/virus_early/virus_config/lib_dark/march"
biaspath = "/Users/gregz/cure/virus_early/virus_config/lib_bias/march"
pixflatpath = "/Users/gregz/cure/virus_early/virus_config/PixelFlats"


# Config dictionary for Amplifier reduction inputs 
config_dict = {'adjust_trace': 'adjust_trace', 'use_trace_ref':'use_trace_ref',
               'use_pixelflat':'use_pixelflat', 'virusconfig':'virusconfig',
               'darkpath':'darkpath','biaspath':'biaspath', 'pixflatpath':'pixflatpath',
               'check_fibermodel':'check_fibermodel', 'check_wave':'check_wave',
               'cosmic_iterations':'cosmic_iterations'}

# Parameter dictionary for Amplifier reduction inputs 
param_dict = { 'fsize':'fs',
               'fibmodel_nbins':'bn', 'wave_nbins':'wbn', 'default_fib':'dfn', 
               'sigma': 'sig', 'power':'pow', 'cont_smooth':'contsmooth',
               'fibmodel_slope':'slope',
               'fibmodel_intercept':'intercept',
               'fibmodel_breakpoint':'breakpoint',
               'fibmodel_step':'fstep','fibmodel_interpkind':'interpk',
               'fiber_group':'fibergroup', 'col_group':'colgroup',
               'wave_res': 'wr'}

param_amp_dict = {'init_lims':'wl', 'specname':'sn', 'dark_mult':'dm',
                  'bias_mult':'bm','collapse_lims':'cwl'}

# Bottom Amplifier for each side 
#VIRUS AND LRS2
Amps = ["LL","RU"]
Sides = ["L", "R"]
Side_dict = {"L": ["LL","LU"], "R": ["RU","RL"]}
Amp_dict = {"LL": ["LU","L"], "RU": ["RL","R"]}

#VIRUS-W Low Res
Amps1 = ["LL"]
Sides1 = ["L"]
Side_dict1 = {"L": ["LL"]}
Amp_dict1 = {"LL": ["LU","L"]}

# VIRUS-W High Res
Amps2 = ["RU"]
Sides2 = ["R"]
Side_dict2 = {"R": ["RU"]}
Amp_dict2 = {"RU": ["RL","R"]}


# Wavelength limits for each side as defined by the bottom amplifier id
virus_wl = {"LL": [3480,5510], "RU": [3480,5510]}
lrs2b_wl = {"LL": [3633,4655], "RU": [4550,7000]}
lrs2r_wl = {"LL": [6425,8440], "RU": [8230,10550]}
virusw_wl = {"LL": [3830,6000],"RU": [4727,5503]}
virusw_lowres_initsol = np.array([155.79487795, -420.19605654,  
                                  2458.0539946 , 3791.40366259])

# Collapsing wavelengths
virus_cwl = {"LL": [4900,5350], "RU": [4900,5350]}
lrs2b_cwl = {"LL": [3633,4655], "RU": [4550,7000]}
lrs2r_cwl = {"LL": [6425,8440], "RU": [8230,10550]}
virusw_cwl = {"LL": [4900,5350], "RU": [4900,5350]}

# Dark multiplier for dark subtraction
virus_dm = {"LL": 1.0, "LU": 1.0, "RU": 1.0, "RL": 1.0}
lrs2b_dm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
lrs2r_dm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
virusw_dm ={"LL": 0.0, "RU": 0.0}

# Bias multiplier for bias subtraction
virus_bm = {"LL": 1.0, "LU": 1.0, "RU": 1.0, "RL": 1.0}
lrs2b_bm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
lrs2r_bm = {"LL": 0.0, "LU": 0.0, "RU": 0.0, "RL": 0.0}
virusw_bm ={"LL": 0.0, "RU": 0.0}

# Bias multiplier for bias subtraction
virus_nf = {"LL": 112, "LU": 112, "RU": 112, "RL": 112}
lrs2b_nf = {"LL": 140, "LU": 140, "RU": 140, "RL": 140}
lrs2r_nf = {"LL": 140, "LU": 140, "RU": 140, "RL": 140}
virusw_nf ={"LL": 267, "RU": 267}

# Name prefix for the normalized spectrum used for the wavelength solution
virus_sn = {"LL": "virus", "RU": "virus"}
lrs2b_sn = {"LL": "lrs2_uv", "RU": "lrs2_orange"}
lrs2r_sn = {"LL": "lrs2_red", "RU": "lrs2_farred"}
virusw_sn = {"LL": "virusw_lowres", "RU": "virusw"}

# Name of the IFUcen file for fiber positions
virus_fn = {"L": ["IFUcen_VIFU",30], "R": ["IFUcen_VIFU",30]}
lrs2b_fn = {"L": ["LRS2_B_UV_mapping.txt",4], "R": ["LRS2_B_OR_mapping.txt",4]}
lrs2r_fn = {"L": ["LRS2_R_NR_mapping.txt",4], "R": ["LRS2_R_FR_mapping.txt",4]}
virusw_fn = {"L": ["virusw_mapping.txt",4],"R": ["virusw_mapping.txt",4]}

# Pixel width in radius over which the fibermodel is defined
virus_fs = 8.
lrs2b_fs = 6.
lrs2r_fs = 6.
virusw_fs = 5.

# Number of bins for fibermodel fit
virus_bn = 31
lrs2b_bn = 31
lrs2r_bn = 31
virusw_bn = 45

# 
virus_slope = 0.000
lrs2b_slope = 0.001
lrs2r_slope = 0.001
virusw_slope = 0.000

# 
virus_intercept = 0.000
lrs2b_intercept = 0.002
lrs2r_intercept = 0.002
virusw_intercept = 0.000

# 
virus_breakpoint = 5.5
lrs2b_breakpoint = 4.
lrs2r_breakpoint = 4.
virusw_breakpoint = 4.

# 
virus_fstep = 4
lrs2b_fstep = 4
lrs2r_fstep = 4
virusw_fstep = 4

# 
virus_interpk = 'linear'
lrs2b_interpk = 'linear'
lrs2r_interpk = 'linear'
virusw_interpk = 'linear'

#
virus_fibergroup = 4
lrs2b_fibergroup = 4
lrs2r_fibergroup = 4
virusw_fibergroup = 4

#
virus_colgroup = 20
lrs2b_colgroup = 40
lrs2r_colgroup = 40
virusw_colgroup = 40

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
virusw_sig = 1.5

# The initial fibermodel power for initializing bins
virus_pow = 2.5
lrs2b_pow = 2.0
lrs2r_pow = 2.0
virusw_pow = 2.8

# Wavelength resolution (for wavelength fit initial guess)
virus_wr = 1.97
lrs2b_wr = 0.5
lrs2r_wr = 1.0
virusw_wr = 0.19

# Dispersion scale for making the Fe/CuFe files
virus_di = {"L": 1.97, "R": 1.97}
lrs2b_di = {"L": 0.5, "R": 1.2}
lrs2r_di = {"L": 1.0, "R": 1.0}
virusw_di = {"L":0.53,"R": 0.19}

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