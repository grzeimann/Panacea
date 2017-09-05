# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 10:39:04 2017

@author: gregz
"""
from args import parse_args
import os.path as op
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import datetime
import warnings
import sys
import traceback
from fiber import Fiber
from fiber_utils import get_indices
from scipy.signal import medfilt
from utils import biweight_location

# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.sans-serif'] = "Meiryo"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "sans-serif"
plt.style.use('seaborn-colorblind')

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

def get_multi_filename(amp):
    return op.join(amp.path, 'multi_%s_%s_%s_%s.fits' 
                             %(amp.specid, amp.ifuslot,
                               amp.ifuid, amp.amp))
def get_fits(amp):
    fn = get_multi_filename(amp)
    try: 
        return fits.open(fn)
    except:
        return None

def write_fits(F):
    try:
        F.writeto(F.filename(),overwrite=True)
    except TypeError:
        F.writeto(F.filename(), clobber=True)
        
def build_sky_list(amps):
    return [amp for amp in amps if (amp.exptime<=900.) and (amp.exptime>1e-8)]

def build_sci_list(amps):
    return [amp for amp in amps if amp.exptime>900.]
    
def reduce_twi(twi, args):
    if op.exists(get_multi_filename(twi)) and args.rerun_twi is False:
        return False
    operations = ['prepare_image', 'get_trace',
                  'get_fibermodel', 'get_wavelength_solution', 
                  'get_fiber_to_fiber']
    for operation in operations:
            execute_function(twi, operation)
    execute_function(twi, 'save_fibmodel')
    image_list = ['image','error']
    spec_list = ['trace', 'wavelength', 'spectrum',  'fiber_to_fiber', 'dead']
    execute_function(twi, 'save', 
                     {'image_list': image_list, 
                      'spec_list': spec_list})
    clean_amp(twi)
    return True

def reduce_sci(sci, args, calpath):
    if op.exists(get_multi_filename(sci)) and args.rerun_sci is False:
        return False
    sci.calpath = calpath
    sci.check_fibermodel = False
    sci.check_wave = False
    operations = ['prepare_image', 'get_trace']
    sci.refit = True
    for operation in operations:
        execute_function(sci, operation)
    sci.refit = False
    execute_function(sci, 'load', {'path':'calpath',
                                   'spec_list':['wavelength', 'spectrum',
                                                'fiber_to_fiber', 'dead']})
    image_list = ['image','error']
    execute_function(sci, 'fiberextract')
    spec_list = ['spectrum', 'wavelength', 'trace', 'fiber_to_fiber']   
    execute_function(sci, 'save', {'image_list': image_list,
                                   'spec_list': spec_list})
    clean_amp(sci)

def get_datetime(amp):
    base = op.basename(amp.filename)
    Y, M, D, H, m, S = [int(s) for s in [base[:4],base[4:6],
                                         base[6:8],base[9:11],
                                         base[11:13],base[13:15]]]
    return datetime.datetime(Y, M, D, H, m, S)

def string_from_datetime(d):
    string = '%04d-%02d-%02d %02d:%02d:%02d' %(d.year, d.month, d.day, d.hour,
                                               d.minute, d.second)
    return string

def get_master_sky(wavelength, spectrum, fiber_to_fiber, path, exptime):
    masterwave = []
    masterspec = []
    for wave, spec, ftf in zip(wavelength,spectrum,fiber_to_fiber):
        masterwave.append(wave)
        masterspec.append(np.where(ftf>1e-8, spec/ftf, 0.0))
    masterwave = np.hstack(masterwave)
    ind = np.argsort(masterwave)
    masterwave[:] = masterwave[ind]
    masterspec = np.hstack(masterspec)
    masterspec[:] = masterspec[ind]
    mastersky = medfilt(masterspec, 281)
    wv = np.arange(masterwave.min(),masterwave.max()+0.05,0.05)
    s = np.zeros((len(wv),2))
    s[:,0] = wv
    s[:,1] = np.interp(wv, masterwave, mastersky / exptime)    
    np.savetxt(op.join(path, 'sky_model.txt'), s)

def reduce_sky(sky, args, trace_list, calpath):
    if op.exists(get_multi_filename(sky)) and args.rerun_sky is False:
        F = get_fits(sky)
        sky.log.info('Getting master sky for %s' %sky.basename)
        return get_master_sky(F['wavelength'].data, F['spectrum'].data, 
                              F['fiber_to_fiber'].data, sky.path, sky.exptime)
    sky.calpath = calpath
    sky.check_fibermodel = False
    sky.check_wave = False
    operations = ['prepare_image']
    for operation in operations:
        execute_function(sky, operation)
        
    get_trace(sky, trace_list)
    execute_function(sky, 'load', {'path':'calpath',
                                   'spec_list':['wavelength', 'spectrum',
                                                'fiber_to_fiber', 'dead']})
    image_list = ['image','error']
    execute_function(sky, 'fiberextract')
    spec_list = ['spectrum', 'wavelength', 'trace', 'fiber_to_fiber']   
    execute_function(sky, 'save', {'image_list': image_list,
                                   'spec_list': spec_list})
    sky.get_master_sky(sky=True)
    wv = np.arange(sky.masterwave.min(),sky.masterwave.max()+0.05,0.05)
    s = np.zeros((len(wv),2))
    s[:,0] = wv
    s[:,1] = np.interp(wv, sky.masterwave, sky.mastersky / sky.exptime)    
    np.savetxt(op.join(sky.path, 'sky_model.txt'), s)
    clean_amp(sky)

def get_interp_weights(amp_in, amp_list):
    timediff = np.zeros((len(amp_list),))
    for i,amp in enumerate(amp_list):
        dt = get_datetime(amp)
        ds = get_datetime(amp_in)
        td = dt - ds
        timediff[i] = td.days * 86400. + td.seconds
    timediff = np.sort(timediff)
    f = np.zeros(timediff.shape)
    w = 1. * f
    for i in np.arange(timediff.shape[0]):
        f[i] = 1.
        w[i] = np.interp(0., timediff, f)
        f[i] = 0.
    return w, np.argsort(timediff)

def get_trace(sky, trace_list):
    sky.log.info('Getting trace from other frames near in time.')
    w, sorted_ind = get_interp_weights(sky, trace_list)
    arr_list = []
    for i in np.arange(len(w)):
        if w[i] > 0.0:
            arr_list.append(w[i] * get_fits(trace_list[sorted_ind[i]])['trace'].data)
    set_trace(sky, np.sum(arr_list, axis=(0,)))


def set_trace(sky, array):
    sky.log.info('Setting trace from other frames near in time.')
    try:
        a = array.shape[0]
    except KeyError:
        sky.log.error('Failed to open extension %s for %s' %('trace', sky.basename))
        return None
    for j in np.arange(a):
        try:
            f = sky.fibers[j]
        except IndexError:    
            f = Fiber(sky.D, j+1, sky.path, sky.filename)                
            sky.fibers.append(f)
        sky.fibers[j].trace = array[j]
    get_indices(sky.image, sky.fibers, sky.fsize)
    
def clean_amp(amp):
    amp.image = None
    amp.back = None
    amp.clean_image = None
    amp.continuum_sub = None
    amp.residual = None
    amp.error = None
    amp.sig = None
    amp.sigwave = None
    amp.error_analysis = None
    amp.fibers = None
    
def subtract_sky_from_sci(sci, sky_list, sky_model_list, wave, use_sci=False):
    sci.log.info('Subtracting sky from %s' %sci.basename)
    if use_sci:
        fn = op.join(sci.path, 'sky_model.txt')
        if True:#not op.exists(fn):
            F = fits.open('panacea/leoI_20131209.fits')
            G = get_fits(sci)
            get_master_sky(G['wavelength'].data, 
                           G['spectrum'].data - F[0].data*sci.exptime,
                           G['fiber_to_fiber'].data, sci.path, sci.exptime)
        wave, sky_model = np.loadtxt(fn, unpack=True)
    else:   
        weight, sorted_ind = get_interp_weights(sci, sky_list)
        arr_list = []
        for i,w in enumerate(weight):
            if w > 0.0:
                arr_list.append(w * sky_model_list[sorted_ind[i]])
        sky_model = np.sum(arr_list, axis=(0,))
    SCI = get_fits(sci)
    sky_subtracted = np.zeros(SCI['spectrum'].data.shape)
    for i, spec in enumerate(SCI['spectrum'].data):
        sky_subtracted[i,:] = (spec - np.interp(SCI['wavelength'].data[i,:],
                                                wave, sky_model*sci.exptime) 
                                      * SCI['fiber_to_fiber'].data[i,:])
    s = fits.ImageHDU(sky_subtracted)
    erase = []
    for i,S in enumerate(SCI):
        if S.header['EXTNAME'] == 'sky_subtracted':
            erase.append(i)
    for i in sorted(erase,reverse=True):
        del SCI[i]
    SCI.append(s)
    SCI[-1].header['EXTNAME'] = 'sky_subtracted'
    write_fits(SCI)            
    return SCI['wavelength'].data, sky_subtracted, wave, sky_model

def make_collapsed_cube(sci, ifucen, ext='sky_subtracted',
                        scale=1.0, seeing=2.0, wlow=5150, 
                        whigh=5250):
    F = get_fits(sci)
    data = np.zeros((F[ext].data.shape[0],))
    for i, v in enumerate(data):
        xl = np.searchsorted(F['wavelength'].data[i,:],wlow,side='left')
        xh = np.searchsorted(F['wavelength'].data[i,:],whigh,side='right')
        data[i] = biweight_location(F[ext].data[i,xl:xh])
    x = np.arange(ifucen[:,0].min()-scale, 
                  ifucen[:,0].max()+scale, scale)
    y = np.arange(ifucen[:,1].min()-scale, 
                  ifucen[:,1].max()+scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)
    d = np.zeros((len(ifucen[:,1]),)+xgrid.shape)
    w = np.zeros((len(ifucen[:,1]),)+xgrid.shape)
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            d[:,j,i]= np.sqrt((ifucen[:,0] - xgrid[j,i])**2 + 
                        (ifucen[:,1] - ygrid[j,i])**2)
            w[:,j,i] = np.exp(-1./2.*(d[:,j,i]/seeing)**2)    
    ws = w.sum(axis=0)
    zgrid = (data[:,np.newaxis,np.newaxis]*w).sum(axis=0)/ws
                
    hdu = fits.PrimaryHDU(np.array(zgrid, dtype='float32'))
    hdu.header['CRPIX1'] = len(x)/2.
    hdu.header['CRPIX2'] = len(y)/2.   
    hdu.header['CRVAL1'] = 0.
    hdu.header['CRVAL2'] = 0.
    hdu.header['CDELT1'] = scale
    hdu.header['CDELT2'] = scale
    erase = []
    for i,S in enumerate(F):
        if S.header['EXTNAME'] == 'collapsed':
            erase.append(i)
    for i in sorted(erase,reverse=True):
        del F[i]
    F.append(hdu)
    F[-1].header['EXTNAME'] = 'collapsed'
    write_fits(F)
        
def main(argv=None):
    args = parse_args(argv)
    ifucen = np.loadtxt(op.join(args.kwargs['virusconfig'], 'IFUcen_files', 
                            args.ifucen_fn['R'][0]), 
                  usecols=[0,1,2], skiprows=args.ifucen_fn['R'][1])
    args.rerun_twi = False
    args.rerun_sci = False
    args.rerun_sky = False
    plot_sky = True
    plot_sci = True
    use_sci = True # first subtract from sky then from sci model
    
    sky_list = build_sky_list(args.sci_list)
    sci_list = build_sci_list(args.sci_list)
    
    for twi in args.twi_list:
        response = reduce_twi(twi, args)
    
    for sci in sci_list:
        reduce_sci(sci, args, args.twi_list[0].path)
    
    trace_list = [sci for sci in sci_list]
    trace_list.insert(0,args.twi_list[0])   
    
    sky_spec_list = []
    if plot_sky:
        wlow = 4800
        whigh = 4900
        fig, ax = plt.subplots(1, figsize=(8,6))
    for sky in sky_list:
        fn = op.join(sky.path, 'sky_model.txt')
        if op.exists(fn):
            w, s = np.loadtxt(fn, unpack=True)
        else:
            reduce_sky(sky, args, trace_list, args.twi_list[0].path)
            w, s = np.loadtxt(op.join(sky.path, 'sky_model.txt'), unpack=True)

        sky_spec_list.append(s)
        if plot_sky:
            xl = np.searchsorted(w, wlow)
            xh = np.searchsorted(w, whigh, side='right')
            ax.plot(w[xl:xh], s[xl:xh], 
                    label='sky:'+string_from_datetime(get_datetime(sky)))
    
    w = np.arange(4700., 5510., 0.05)                
    if plot_sci:
        wlow = 4800
        whigh = 4900
        fig1, ax1 = plt.subplots(1, figsize=(8,6))
        fibers = [76]
    for sci in sci_list:        
        wave, skys, wvs, skysp = subtract_sky_from_sci(sci, sky_list, sky_spec_list, w,
                                           use_sci=use_sci)
        make_collapsed_cube(sci, ifucen[:,1:3][::-1,:])
        if plot_sci:
            for fiber in fibers:
                xl = np.searchsorted(wave[fiber,:], wlow)
                xh = np.searchsorted(wave[fiber,:], whigh, side='right')
                fstr = '%03d:' %fiber
                ax1.plot(wave[fiber,xl:xh], skys[fiber,xl:xh]/sci.exptime, 
                        label=fstr+string_from_datetime(get_datetime(sci)))
        if plot_sky:
            xl = np.searchsorted(wvs, wlow)
            xh = np.searchsorted(wvs, whigh, side='right')
            ax.plot(wvs[xl:xh], skysp[xl:xh], 
                    label='sci:'+string_from_datetime(get_datetime(sci)))
    if plot_sci:
        ax1.set_xlim([wlow,whigh])
        ax1.set_ylim([0.00, 0.1])
        ax1.legend(loc='best', fancybox=True, framealpha=0.5)
        ax1.tick_params(labelsize=10)
        ax1.set_ylabel(r'Sci Brightness (e-/s)', fontsize=12)
        ax1.set_xlabel(r'Wavelength ($\AA$)', fontsize=12)
        fig1.savefig('sci_spec_%s.png' %args.twi_list[0].date.date(), dpi=150)
    if plot_sky:
        ax.set_xlim([wlow,whigh])
        ax.set_ylim([0.00, 0.05])
        ax.legend(loc='best', fancybox=True, framealpha=0.5)
        ax.tick_params(labelsize=10)
        ax.set_ylabel(r'Sky Brightness (e-/s)', fontsize=12)
        ax.set_xlabel(r'Wavelength ($\AA$)', fontsize=12)
        fig.savefig('sky_spec_%s.png' %args.twi_list[0].date.date(), dpi=150)

    
    
            
if __name__ == '__main__':
    main()    
    
    
            
    