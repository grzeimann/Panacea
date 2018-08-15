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
try:
    from pyhetdex.cure.distortion import Distortion
    pyhetdex_flag = True
except:
    print('Pyhetdex not installed.  Continuing on.')
    pyhetdex_flag = False
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

def get_ifucenfile(args, side, ifuid, amp=None):
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
            if (ifuid == '004'):
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
                if ifuid in ['003','005','008']:
                    ifucen[224:,:] = ifucen[-1:223:-1,:]
                if ifuid == '007':
                    ifucen[37,1:3], ifucen[38,1:3] = ifucen[38,1:3], ifucen[37,1:3]
                if ifuid == '025':
                    ifucen[208,1:3], ifucen[213,1:3] = ifucen[213,1:3], ifucen[208,1:3]
                if ifuid == '030':
                    ifucen[445,1:3], ifucen[446,1:3] = ifucen[446,1:3], ifucen[445,1:3]
                if ifuid == '038':
                    ifucen[302,1:3], ifucen[303,1:3] = ifucen[303,1:3], ifucen[302,1:3]
                if ifuid == '041':
                    ifucen[251,1:3], ifucen[252,1:3] = ifucen[252,1:3], ifucen[251,1:3]
    else:
        ifucen = np.loadtxt(op.join(args.kwargs['virusconfig'], 'IFUcen_files', 
                            args.ifucen_fn[side][0]), 
                  usecols=[0,1,2], skiprows=args.ifucen_fn[side][1])   
                  
    if amp is None:
        return ifucen
    else:
        if args.instr=="virus" and args.kwargs['use_trace_ref']:
            if amp=="LL":
                return ifucen, ifucen[112:224,1:3][::-1,:]
            if amp=="LU":
                return ifucen, ifucen[:112,1:3][::-1,:]
            if amp=="RL":
                return ifucen, ifucen[224:336,1:3][::-1,:]               
            if amp=="RU":
                return ifucen, ifucen[336:,1:3][::-1,:] 
        if args.instr=='virusw':
            return ifucen, ifucen[:,1:3][::-1,:]
        if args.instr=='lrs2':
            if amp=="LL":
                return ifucen, ifucen[140:,1:3][::-1,:]
            if amp=="LU":
                return ifucen, ifucen[:140,1:3][::-1,:]
            if amp=="RL":
                return ifucen, ifucen[:140,1:3][::-1,:]               
            if amp=="RU":
                return ifucen, ifucen[140:,1:3][::-1,:] 
                
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
            ifucen, temp = get_ifucenfile(args, amp.amp[0], 
                                          amp.ifuid, amp.amp)
            amp.ifupos = temp
            image_list.append('ifupos')
            execute_function(amp, 'save', {'image_list':image_list,
                                           'spec_list':['trace','wavelength',
                                                        'spectrum',
                                                        'fiber_to_fiber',
                                                        'dead']})
    if args.reduce_sci:
        if args.twipath is not None:
            calpath = args.twipath
        else:
            calpath = args.twi_list[0].path
        reduce_list = ['sci_list']
        for _list in reduce_list:
            if hasattr(args, _list):
                for amp in getattr(args, _list):
                    if args.use_other_sky and _list=='sci_list':
                        try:
                            amp.skypath = args.sky_list[0].path
                        except:
                            amp.log.error('Could not set skypath, using sky in frame instead.')
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
                    if amp.adjust_trace or args.trace_from_sci:
                        amp.refit=True
                        execute_function(amp, 'get_trace')
                        amp.refit=False
                        
                        
                    if args.sci_operations['remeasure_fibermodel']:  
                        amp.refit=True
                        amp.check_fibermodel=True                                          
                        execute_function(amp, 'get_fibermodel')
                        execute_function(amp, 'save_fibmodel')
                        amp.refit=False
                    else:
                        execute_function(amp, 'get_fibermodel')
                    if args.refit_fiber_to_fiber:
                        amp.refit=True
                        execute_function(amp, 'get_fiber_to_fiber')
                        amp.refit=False
                    if args.adjust_ftf:
                        execute_function(amp, 'get_fiber_to_fiber')
                    execute_function(amp, 'fiberextract')
                    if args.sci_operations['sky_subtraction']:                                            
                        execute_function(amp, 'sky_subtraction')
                        image_list.append('clean_image')
                        image_list.append('flat_image')
                        image_list.append('continuum_sub')
                        image_list.append('residual')
                        spec_list.append('sky_subtracted')
                        spec_list.append('sky_spectrum')
                    if amp.cosmic_iterations>0.0:
                        execute_function(amp, 'clean_cosmics')
                        execute_function(amp, 'fiberextract')
                        if args.sci_operations['sky_subtraction']:                                            
                            execute_function(amp, 'sky_subtraction')
                    if args.instr in ['virus', 'virusw', 'lrs2']:
                        ifucen, temp = get_ifucenfile(args, amp.amp[0], 
                                            amp.ifuid, amp.amp)
                        amp.ifupos = temp
                        image_list.append('ifupos')
                    if args.sci_operations['significance_map']:
                        execute_function(amp, 'get_significance_map')
                        image_list.append('sig')
                        image_list.append('sigwave')
                    if args.kwargs['make_model_image']:
                        image_list.append('fibmodel_image')
                    if args.sci_operations['sky_subtraction']:                                            
                        execute_function(amp, 'make_error_analysis')
                        image_list.append('error_analysis')
                    execute_function(amp, 'save', {'image_list':image_list,
                                                   'spec_list':spec_list})
                    amp.image = None
                    amp.back = None
                    amp.clean_image = None
                    amp.continuum_sub = None
                    amp.residual = None
                    amp.error = None
                    #amp.sig = None
                    #amp.sigwave = None
                    #amp.error_analysis = None
                    #amp.fibers = None
            
            
    if args.combine_reductions:
        paths = np.array([amp.path for amp in args.sci_list])
        unique_paths = np.unique([amp.path for amp in args.sci_list])
        if pyhetdex_flag:
            D = get_distortion_file(args)
        for up in unique_paths:
            loc = np.where(up==paths)[0][0]
            print(loc)
            spec = Spectrograph(up, args.sci_list[loc].specid, 
                                args.sci_list[loc].ifuslot, 
                                args.sci_list[loc].ifuid,
                                args.sci_list[loc].basename, 
                                N=args.sci_list[loc].N, D=args.sci_list[loc].D,
                                nfib=len(args.sci_list[loc].fibers),
                                scale=args.scale,
                                side_dict = args.side_dict,
                                sides = args.sides, header=args.sci_list[loc].header)
                    
            for side in args.side_dict:
                ids =  [s for s in np.where(up==paths)[0]
                          if args.sci_list[s].amp[0]==side]
                spec.wavelim = args.sci_list[ids[0]].init_lims
                spec.collapselim = args.sci_list[ids[0]].collapse_lims
                spec.disp = args.disp[side]
                if not args.limited_output_files:
                    if args.kwargs['make_model_image']:
                        execute_function(spec, 'write_spectrograph_image',
                                     {  'side':side, 'ext':'fibmodel_image', 
                                     'prefix':'model_'})
                    
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
                if pyhetdex_flag:
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
                if not args.limited_output_files:
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
                if args.sci_operations['sky_subtraction']:                                            
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