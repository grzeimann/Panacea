# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 11:33:22 2016

@author: gregz
"""

from amplifier import Amplifier
import time

#"075","076","083","085","086","093","094","095","096","103","104","105","106"
ifuslots = ["106"]
amps = ["LL"]
base = "/Users/gregz/cure/virus_raw/20161206/virus/virus0000001/exp02/virus/20161206T001914.7_"

for ifu in ifuslots:
    for amp in amps:
        t1 = time.time()
        name = "%s%s_twi.fits" %(ifu, amp)
        full_name = base + name
        twi = Amplifier(full_name, '/Users/gregz/cure/virus_reductions/20161206/twi',
                  calpath='/Users/gregz/cure/virus_reductions/20161206/twi', 
                  debug=True, refit=True, dark_mult=0.0)
        #twi.load_fibers()
        twi.get_fiber_to_fiber(use_default_profile=False, 
                               init_lims=[3490,5500], interactive=True,
                               check_fibermodel=True, check_wave=True)
        twi.save_fibers()
        t2  = time.time()
        print("Time Taken for %s_%s: %0.2f s" %(ifu, amp, t2-t1))
