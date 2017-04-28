# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:02:22 2017

@author: gregz
"""

import matplotlib
matplotlib.use('agg')
from args import parse_args
import time
import sys
import traceback
import numpy as np
from spectrograph import Spectrograph
import os.path as op
from pyhetdex.cure.distortion import Distortion
import warnings


def execute_function(obj, call, kwargs={}):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            func = getattr(obj, call)
            func(**kwargs)
    except:
        obj.log.error('Error occured while running %s on %s' %(call, obj.basename))
        obj.log.error(sys.exc_info()[0])
        obj.log.error(traceback.print_exc(file=sys.stdout))

def get_ifucenfile(args, side, ifuid):
    if args.instr == "virus":
        if not args.kwargs['use_trace_ref']:
            ifucen = np.loadtxt(op.join(args.kwargs['virusconfig'], 
                                        'IFUcen_files', 
                                        args.ifucen_fn[side][0]
                                        + ifuid
                                        + '.txt'), 
                                        usecols=[0,1,2,4], 
                                   skiprows=args.ifucen_fn[side][1])
            
        else:
            if ifuid == '004':
                ifucen = np.loadtxt(op.join(args.kwargs['virusconfig'],
                                        'IFUcen_files',
                                        'IFUcen_HETDEX_reverse_R.txt'),
                                        usecols=[0,1,2,4],
                                   skiprows=args.ifucen_fn[side][1])
                ifucen[224:,:] = ifucen[-1:223:-1,:]
            else:
                ifucen = np.loadtxt(op.join(args.kwargs['virusconfig'],
                                        'IFUcen_files',
                                        'IFUcen_HETDEX.txt'),
                                        usecols=[0,1,2,4],
                                   skiprows=args.ifucen_fn[side][1])
    else:
        ifucen = np.loadtxt(op.join(args.kwargs['virusconfig'], 'IFUcen_files', 
                            args.ifucen_fn[side][0]), 
                  usecols=[0,1,2], skiprows=args.ifucen_fn[side][1])   
    return ifucen
    
def get_distortion_file(args):
    return Distortion(op.join(args.kwargs['virusconfig'], 'DeformerDefaults', 
                                        'mastertrace_twi_027_L.dist')) 

def main():
    t1 = time.time()
    args = parse_args()
    if args.reduce_twi:
        for amp in args.twi_list:
            operations = ['prepare_image',  'subtract_background', 'get_trace',
                          'get_fibermodel', 'get_wavelength_solution', 
                          'get_fiber_to_fiber']
            for operation in operations:
                if args.twi_operations[operation]:
                    execute_function(amp, operation)
            execute_function(amp, 'save_fibmodel')
            image_list = ['image','error']
            if args.twi_operations['subtract_background']:
                image_list.insert(1, 'back')
            execute_function(amp, 'save', {'image_list':image_list,
                                           'spec_list':['trace','wavelength',
                                                        'spectrum',
                                                        'fiber_to_fiber',
                                                        'dead']})
    if args.reduce_sci:
        calpath = args.twi_list[0].path
        for amp in args.sci_list:
            amp.calpath = calpath
            amp.check_fibermodel=False
            amp.check_wave=False
            execute_function(amp, 'prepare_image')
            execute_function(amp, 'load', {'path':'calpath',
                                           'spec_list':['trace','wavelength',
                                                        'spectrum',
                                                        'fiber_to_fiber',
                                                        'dead']})
            image_list = ['image','error']
            spec_list = ['spectrum','wavelength','trace','fiber_to_fiber',
                         'twi_spectrum']                                           
            if args.sci_operations['subtract_background']:                                            
                execute_function(amp, 'subtract_background')
                image_list.append('back')
            if amp.adjust_trace:
                amp.refit=True
                execute_function(amp, 'get_trace')
                amp.refit=False
            execute_function(amp, 'load_fibmodel', {'path':'calpath'})
            if args.refit_fiber_to_fiber:
                amp.refit=True
                execute_function(amp, 'get_fiber_to_fiber')
                amp.refit=False
            if args.sci_operations['sky_subtraction']:                                            
                execute_function(amp, 'sky_subtraction')
                image_list.append('clean_image')
                image_list.append('continuum_sub')
                image_list.append('residual')
                spec_list.append('sky_subtracted')
                spec_list.append('corrected_spectrum')
                spec_list.append('corrected_sky_subtracted')
            if amp.cosmic_iterations>0.0:
                execute_function(amp, 'clean_cosmics')
                execute_function(amp, 'fiberextract')
                if args.sci_operations['sky_subtraction']:                                            
                    execute_function(amp, 'sky_subtraction')
            execute_function(amp, 'save', {'image_list':image_list,
                                           'spec_list':spec_list})
            amp.image = None
            amp.back = None
            amp.clean_image = None
            amp.continuum_sub = None
            amp.residual = None
            amp.error = None
            
            
    if args.combine_reductions:
        paths = np.array([amp.path for amp in args.sci_list])
        unique_paths = np.unique([amp.path for amp in args.sci_list])
        D = get_distortion_file(args)
        for up in unique_paths:
            loc = np.where(up==paths)[0][0]
            spec = Spectrograph(up, args.sci_list[loc].specid, 
                                args.sci_list[loc].ifuslot, 
                                args.sci_list[loc].ifuid,
                                args.sci_list[loc].basename, 
                                N=args.sci_list[loc].N, D=args.sci_list[loc].D,
                                nfib=len(args.sci_list[loc].fibers),
                                scale=args.scale)
                    
            for side in spec.side_dict:
                ids =  [s for s in np.where(up==paths)[0]
                          if args.sci_list[s].amp[0]==side]
                spec.wavelim = args.sci_list[ids[0]].init_lims
                spec.collapselim = args.sci_list[ids[0]].collapse_lims
                spec.disp = args.disp[side]
                execute_function(spec, 'write_spectrograph_image',
                                 {'side':side, 'ext':'image', 'prefix':''})
                execute_function(spec, 'write_spectrograph_image',
                                 {'side':side, 'ext':'error', 'prefix':'e.'})
                if args.sci_operations['sky_subtraction']:                                            

                    execute_function(spec, 'write_spectrograph_image',
                                     {'side':side, 'ext':'clean_image', 
                                      'prefix':'S'})
                    execute_function(spec, 'write_spectrograph_image',
                                     {'side':side, 'ext':'error', 
                                      'prefix':'e.S'})
                execute_function(spec, 'write_new_distfile', {'D':D,
                                                              'side':side})
             
                if args.instr is not "virus":
                    ifucen = get_ifucenfile(args, side, 
                                            args.sci_list[loc].ifuid)
                    spec.ifucen = ifucen
                    spec.scale = args.scale
                    execute_function(spec, 'write_fiberextract',
                                 {'side':side, 'ext':'spectrum', 
                                 'prefix':'Fe'})
                    execute_function(spec, 'write_cube', {'ext':'spectrum',
                                                          'prefix':['CuFe','CoFe'],
                                                          'side':side})
                    if args.sci_operations['sky_subtraction']:                                            
                        execute_function(spec, 'write_fiberextract',
                                         {'side':side, 'ext':'sky_subtracted', 
                                          'prefix':'FeS'})
                        execute_function(spec, 'write_cube', {'ext':'sky_subtracted',
                                                              'prefix':['CuFeS','CoFeS'],
                                                              'side':side})
            if args.instr is "virus":
                ifucen = get_ifucenfile(args, side, 
                                            args.sci_list[loc].ifuid)
                spec.ifucen = ifucen
                spec.scale = args.scale
                for side in spec.side_dict:
                    execute_function(spec, 'write_fiberextract',
                                 {'side':side, 'ext':'spectrum', 
                                 'prefix':'Fe'})
                execute_function(spec, 'write_cube', {'ext':'spectrum',
                                                      'prefix':['CuFe','CoFe']})
                for side in spec.side_dict:
                    execute_function(spec, 'write_fiberextract',
                                 {'side':side, 'ext':'twi_spectrum', 
                                 'prefix':'FeT'})
                execute_function(spec, 'write_cube', {'ext':'twi_spectrum',
                                                      'prefix':['CuFeT','CoFeT']})
                for side in spec.side_dict:
                    execute_function(spec, 'write_fiberextract',
                                 {'side':side, 'ext':'sky_subtracted', 
                                 'prefix':'FeS'})
                execute_function(spec, 'write_cube', {'ext':'sky_subtracted',
                                                      'prefix':['CuFeS','CoFeS']})            
       
    t2 = time.time()
    print("Time Taken: %0.3f" %(t2-t1))

if __name__ == '__main__':
    main() 