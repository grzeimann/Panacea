# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:02:22 2017

@author: gregz
"""

from args import parse_args
import time
import os.path as op
import sys
import traceback

def execute_function(amp, call, kwargs={}):
    try:
        func = getattr(amp, call)
        func(**kwargs)
    except:
        amp.log.error('Error occured while running %s on %s' %(call, amp.basename))
        amp.log.error(sys.exc_info()[0])
        amp.log.error(traceback.print_exc(file=sys.stdout))

        

def main():
    t1 = time.time()
    args = parse_args()
    if args.reduce_twi:
        for amp in args.twi_list:
            execute_function(amp, 'get_fiber_to_fiber')
            execute_function(amp, 'save_fibers', {'cal':True, 'fibmodel':True})
    if args.reduce_sci:
        calpath = args.twi_list[0].path
        for amp in args.sci_list:
            amp.calpath = calpath
            amp.check_fibermodel=False
            amp.check_wave=False
            execute_function(amp, 'prepare_image')
            execute_function(amp, 'load_fibers', {'path':'calpath',
                                                  'cal':True, 'fibmodel':True})
            if amp.adjust_trace:
                amp.refit=True
                execute_function(amp, 'get_trace')
                amp.refit=False
            if args.refit_fiber_to_fiber:
                amp.refit=True
                execute_function(amp, 'get_fiber_to_fiber')
                amp.refit=False
            execute_function(amp, 'sky_subtraction')
            execute_function(amp, 'save', {'image_list':['image','clean_image', 
                                                         'continuum_sub', 
                                                         'residual','error']})
    t2 = time.time()
    print("Time Taken: %0.3f" %(t2-t1))

if __name__ == '__main__':
    main() 