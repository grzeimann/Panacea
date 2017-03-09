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
from amplifier import Amplifier
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
    for spec in args.specid:
        spec_ind_twi = np.where(args.twi_df['Specid'] == spec)[0]
        for amp in config.Amps:
            amp_ind_twi = np.where(args.twi_df['Amp'] == amp)[0]
            twi_sel = np.intersect1d(spec_ind_twi, amp_ind_twi)
            for ind in twi_sel:
                if args.debug:
                    print("Working on Cal for %s, %s" %(spec, amp))                    
                twi1 = Amplifier(args.twi_df['Files'][ind],
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=True, refit=True, dark_mult=0.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir, 
                                 specname=args.specname[amp],
                                 use_pixelflat=(args.pixelflats<1),
                                 init_lims=args.wvl_dict[amp], 
                                 check_fibermodel=True, check_wave=True,
                                 fsize=args.fsize, 
                                 fibmodel_nbins=args.fibmodel_bins,
                                 sigma=args.fibmodel_sig,
                                 power=args.fibmodel_pow)
                twi1.get_trace()
                fname = op.join(args.configdir,'Fiber_Locations',
                                args.twidir_date[0],
                                'fiber_loc_%s_%s_%s_%s.txt' %(
                                args.twi_df['Specid'][ind],
                                args.twi_df['Ifuslot'][ind],
                                args.twi_df['Ifuid'][ind],
                                amp))
                if len(twi1.fibers) == args.nfibers[amp]:
                    print("The correct number of fibers, %i, were found." 
                          %args.nfibers[amp])
                    Y = np.zeros((len(twi1.fibers),2))
                    col = int(twi1.D*config.frac)
                    for i,fiber in enumerate(twi1.fibers):
                        Y[i,0] = fiber.trace[col]
                    print("Printing file for fibers.")
                    create_fiber_file(fname, config.frac, Y)
                else:
                    print("Only found %i fibers instead of the expected %i"
                          %(len(twi1.fibers), args.nfibers[amp]))
                    Y = find_missing_fibers(config.frac, twi1)
                    create_fiber_file(fname, config.frac, Y)
                twi2 = Amplifier(args.twi_df['Files'][ind].replace(amp, 
                                                      config.Amp_dict[amp][0]),
                                 args.twi_df['Output'][ind],
                                 calpath=args.twi_df['Output'][ind], 
                                 debug=True, refit=True, dark_mult=0.0,
                                 darkpath=args.darkdir, biaspath=args.biasdir,
                                 virusconfig=args.configdir, 
                                 specname=args.specname[amp],
                                 use_pixelflat=(args.pixelflats<1),
                                 init_lims=args.wvl_dict[amp], 
                                 check_fibermodel=True, check_wave=True,
                                 fsize=args.fsize, 
                                 fibmodel_nbins=args.fibmodel_bins,
                                 sigma=args.fibmodel_sig,
                                 power=args.fibmodel_pow)
                twi2.get_trace()
                folder = op.join(args.configdir,'Fiber_Locations',
                                args.twidir_date[0])
                mkpath(folder)
                fname = op.join(args.configdir,'Fiber_Locations',
                                args.twidir_date[0],
                                'fiber_loc_%s_%s_%s_%s.txt' %(
                                args.twi_df['Specid'][ind],
                                args.twi_df['Ifuslot'][ind],
                                args.twi_df['Ifuid'][ind],
                                config.Amp_dict[amp][0]))
                               
                if len(twi2.fibers) == args.nfibers[config.Amp_dict[amp][0]]:
                    print("The correct number of fibers, %i, were found." 
                          %args.nfibers[config.Amp_dict[amp][0]])
                    Y = np.zeros((len(twi2.fibers),2))
                    col = int(twi2.D*config.frac)
                    for i,fiber in enumerate(twi2.fibers):
                        Y[i,0] = fiber.trace[col]
                    print("Printing file for fibers.")
                    create_fiber_file(fname, config.frac, Y)
                else:
                    print("Only found %i fibers instead of the expected %i"
                          %(len(twi2.fibers), 
                            args.nfibers[config.Amp_dict[amp][0]]))
                    Y = find_missing_fibers(config.frac, twi2)
                    create_fiber_file(fname, config.frac, Y)
                    
if __name__ == '__main__':
    main()    