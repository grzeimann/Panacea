#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Defining default fiber location 
------------------
Built for the VIRUS instrument as well as LRS2 on HET

Incomplete Documentation

"""

import pyds9 
import config
import numpy as np
import os.path as op
from args import parse_args
from distutils.dir_util import mkpath


def create_fiber_file(fname, frac, array):
    np.savetxt(fname, array, fmt = "%1.3f %i", 
               header='%0.3f'%frac)

def find_missing_fibers(frac, amp_obj):
    col = int(amp_obj.D * frac)
    init_fibers = ['%0.3f' %(fiber.trace[col]+1.) for fiber in amp_obj.fibers]
    ds9 = pyds9.DS9()
    ds9.set_np2arr(amp_obj.image)
    for i,fiber in enumerate(amp_obj.fibers):
        ds9.set('regions command "circle %0.3f %s 2 # color=red"'
                %(col, init_fibers[i]))
    ds9.set('regions format xy')
    raw_input('Hit ENTER when done marking missing fibers...')
    strds9 = ds9.get('regions list')
    str_fib = strds9.split('\n')
    Y = np.zeros((len(str_fib),2))
    for i,strf in enumerate(str_fib):
        s = '%0.3f' %float(strf.split()[1])
        Y[i,0] = float(s)-1.
        if s not in init_fibers:
            Y[i,1] = 1
    Y = Y[Y[:,0].argsort(),:]
    return Y
    
def main():
    args = parse_args()
    for twi in args.twi_list:
        twi.log.info("Working on Cal for %s, %s" %(twi.ifuslot, twi.amp)) 
        twi.use_trace_ref = False                   
        twi.get_trace()
        folder = op.join(args.kwargs['virusconfig'],'Fiber_Locations',
                                '%04d%02d%02d' 
                                %(twi.date.year,twi.date.month, twi.date.day))
        mkpath(folder)              
        fname = op.join(folder, 'fiber_loc_%s_%s_%s_%s.txt' %(
                                twi.specid, twi.ifuslot, twi.ifuid,
                                twi.amp))
    
        if len(twi.fibers) == args.nfibers[twi.amp]:
            twi.log.info("The correct number of fibers, %i, were found." 
                          %args.nfibers[twi.amp])
            Y = np.zeros((len(twi.fibers),2))
            col = int(twi.D*config.frac)
            for i,fiber in enumerate(twi.fibers):
                Y[i,0] = fiber.trace[col]
            twi.log.info("Printing file for fibers.")
            create_fiber_file(fname, config.frac, Y)
        else:
            twi.log.info("Only found %i fibers instead of the expected %i"
                          %(len(twi.fibers), args.nfibers[twi.amp]))
            Y = find_missing_fibers(config.frac, twi)
            create_fiber_file(fname, config.frac, Y)
                    
if __name__ == '__main__':
    main()    