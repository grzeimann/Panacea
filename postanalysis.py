# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:55:42 2017

@author: gregz
"""

import argparse as ap
import glob
import os.path as op
from astropy.io import fits
import matplotlib.pyplot as plt
from pyhetdex.het.fplane import FPlane
from pyhetdex.coordinates.tangent_projection import TangentPlane as TP
import numpy as np
from utils import biweight_location, is_outlier, biweight_midvariance
from photutils import Background2D, SigmaClip, BiweightLocationBackground
from photutils import detect_sources, deblend_sources, DAOStarFinder
from photutils import EllipticalAperture
from astropy.visualization import AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from scipy.signal import medfilt2d
from astroquery.skyview import SkyView
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy.nddata import Cutout2D
from astroquery.vizier import Vizier 
import astropy.units as unit 
import astropy.coordinates as coord
import astropy.units as aunit
import logging
import sys
from astrometry import Astrometry
from fiber import Fiber
from fiber_utils import get_indices
from collections import Counter
from scipy.cluster.vq import kmeans2
from astropy.convolution import convolve, Gaussian1DKernel

sys_rot = 1.20

# common parts of the argument parser
astro_parent = ap.ArgumentParser(add_help=False)
_group_astro = astro_parent.add_mutually_exclusive_group()
_group_astro.add_argument("--ifuslots", nargs='?', type=str, 
                            help='''IFUSLOT for processing. 
                            Ex: "075,076".''', default = None)
_group_astro.add_argument("--ifuslot_file", nargs='?', type=str, 
                            help='''IFUSLOT file. 
                            
                            Ex: "good_ifuslots.txt".''', default = None)
astro_parent.add_argument("--folder", nargs='?', type=str, 
                            help='''Reduction Folder
                            Ex: "reductions".''', default = "reductions")
astro_parent.add_argument("--instr", nargs='?', type=str, 
                            help='''Instrument to process. 
                            Default: "virus"
                            Ex: "camra" for lab data,
                                "lrs2" for lrs2.''', default = "virus")
astro_parent.add_argument("-sd","--scidir_date", nargs='?', type=str,
                            help='''Science Directory Date.    
                            Ex: \"20160412\"''', default=None)
astro_parent.add_argument("-so","--scidir_obsid", nargs='?', type=str,
                            help='''Science Directory ObsID.   
                            Ex: \"3\" or \"102\"''', default=None)                            
astro_parent.add_argument("-se","--scidir_expnum", nargs='?', type=str,
                            help='''Science Directory exposure number.
                            Ex: \"1\" or \"05\"''', default=None)              
astro_parent.add_argument("--fplane", nargs='?', type=str, 
                            help='''Fplane file
                            Ex: "fplane.txt".''', default = 'fplane.txt')


def setup_logging():
    '''Set up a logger for analysis with a name ``shot``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    log = logging.getLogger('shot')
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        level = logging.INFO
       
        fmt = logging.Formatter(fmt)
    
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)
    
        log = logging.getLogger('shot')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log

    
def get_shot_info_from_input(args):
    '''
    Command line or argument input includes scidir_date, scidir_obsid, and
    scidir_expnum.  This input then points to different shots.  This function
    returns the list of dates, obsid and exposure numbers for that input.
    '''    
    log = setup_logging()

    obs = 'sci'                         
    labels = ['dir_date', 'dir_obsid', 'dir_expnum']

    for label in labels[:2]:
        if getattr(args, obs+label) is None:
            log.error('%s%s was not provided' %(obs, label))
            sys.exit(1) 
        if not isinstance(getattr(args, obs+label), list):
            setattr(args, obs+label, 
                    getattr(args, obs+label).replace(" ", "").split(','))
    if getattr(args, obs+labels[2]) is not None \
       and not isinstance(getattr(args, obs+labels[2]), list):
        setattr(args, obs+labels[2], 
                getattr(args, obs+labels[2]).replace(" ", "").split(','))
    shotinfo_list = []
    for date in getattr(args, obs+labels[0]):
        for obsid in getattr(args, obs+labels[1]):
            if getattr(args, obs+labels[2]) is not None:  
                exposures = getattr(args, obs+labels[2])
            else:
                folder = op.join(date, args.instr,
                                 "{:s}{:07d}".format(args.instr, 
                                                     int(obsid)))
                exposures = [op.basename(exp)[-2:] for exp in 
                                   glob.glob(op.join(args.folder, folder,'*'))] 
            for expnum in exposures:
                shotinfo_list.append([date, obsid, expnum])
    return shotinfo_list 

def get_astrometry_info_from_input(args, ifuslot):
    '''
    Command line or argument input includes scidir_date, scidir_obsid, and
    scidir_expnum.  This input then points to different shots.  This function
    returns the list of dates, obsid and exposure numbers for that input.
    '''    
    log = setup_logging()
    
    astrom_info = ['RA', 'DEC', 'PA', 'RHO', 'EXPTIME','OSCANMN','OSCANSTD']
    obs = 'sci'                         
    labels = ['dir_date', 'dir_obsid', 'dir_expnum']

    for label in labels[:2]:
        if getattr(args, obs+label) is None:
            log.error('%s%s was not provided' %(obs, label))
            sys.exit(1) 
        if not isinstance(getattr(args, obs+label), list):
            setattr(args, obs+label, 
                    getattr(args, obs+label).replace(" ", "").split(','))
    if getattr(args, obs+labels[2]) is not None \
       and not isinstance(getattr(args, obs+labels[2]), list):
        setattr(args, obs+labels[2], 
                getattr(args, obs+labels[2]).replace(" ", "").split(','))
    astrominfo_list = []
    for date in getattr(args, obs+labels[0]):
        for obsid in getattr(args, obs+labels[1]):
            if getattr(args, obs+labels[2]) is not None:  
                exposures = getattr(args, obs+labels[2])
            else:
                folder = op.join(date, args.instr,
                                 "{:s}{:07d}".format(args.instr, 
                                                     int(obsid)))
                exposures = [op.basename(exp)[-2:] for exp in 
                                   glob.glob(op.join(args.folder, folder,'*'))] 
            for expnum in exposures:
                path = build_path(args.folder, args.instr, date, obsid, expnum)
                fn = op.join(path, 'multi_*_%s_*_LL.fits' %ifuslot)
                fn = glob.glob(fn)
                if not fn:
                    log.error('No file found at %s for %s' %(path, ifuslot))
                    sys.exit(1)
                F = fits.open(fn[0])
                l = []
                for astrom in astrom_info:
                    l.append(F[0].header[astrom])
                astrominfo_list.append(l)
    return astrominfo_list 


def build_path(reduction_folder, instr, date, obsid, expn):
    folder = op.join(date, instr, "{:s}{:07d}".format(instr, int(obsid)), 
                     "exp{:02d}".format(int(expn)), instr)
    return op.join(reduction_folder, folder)
                                     
 
def read_in_multifiles(reduction_folder, date, obsid, expn, ifuslot, instr, 
                       amps=None, extensions=None):
    '''
    Read in the multi_* files from panacea reduction
    
    Return arrays
    '''
    log = setup_logging()
    if extensions is None:
        log.error('Extensions were not provided for input')
        sys.exit(1)   
        
    list_outer = []
    for i, ext in enumerate(extensions):
        list_outer.append([])
    list_inner = []
    for i, ext in enumerate(extensions):
        list_inner.append([])
    folder = op.join(date, instr, "{:s}{:07d}".format(instr, int(obsid)), 
                     "exp{:02d}".format(int(expn)), instr)
    if amps is None:
        files = sorted(glob.glob(op.join(reduction_folder, folder, 
                                     'multi_*_%s*' %ifuslot)))
    else:
        files = []
        for amp in amps:
            files.append(glob.glob(op.join(reduction_folder, folder, 
                                     'multi_*_%s_*%s.fits' %(ifuslot,amp)))[0])
    if not files:
        log.error('No Files were found for ifuslot %s on %s, %s, %s' 
                  %(ifuslot, date, obsid, expn))
        sys.exit(1)
    for fn in files:
        F = fits.open(fn)
        for i, ext in enumerate(extensions):
            list_inner[i].append(F[ext].data*1.)
            del F[ext].data
        F.close()
    for i, ext in enumerate(extensions):
        list_outer[i].append(np.vstack(list_inner[i]))
    return list_outer

def fiber_to_fiber_shot(shot, folder, instr, ifuslots, amps_list,
                        wave_low=3470., wave_high=5530., 
                        wave_step=1., smooth_scale=75,
                        func=medfilt2d, write_to_file=False):
    '''
    Parameters
    ----------
    args : list of strings, optional
        command line
    '''
    log = setup_logging()

    log.info('Calculating Fiber to Fiber for date/obsid/expn: %s, %s, %s'
             %(shot[0], shot[1], shot[2]))

    wave_list,twi_list = [], []
    for k,ifuslot in enumerate(ifuslots):
        wave, twi = read_in_multifiles(folder, shot[0], shot[1], 
                                       shot[2], ifuslot, instr, amps=amps_list[k], 
                                       extensions=['wavelength',
                                                   'twi_spectrum'])
        wave_list.append(wave)
        twi_list.append(twi)
    std_wave = np.arange(wave_low,wave_high+wave_step,wave_step)
    avg_spec_list = []
    for wave, twi in zip(wave_list, twi_list):
        wave, twi = wave[0], twi[0]
        a,b = wave.shape
        S = np.zeros((a,len(std_wave)))
        for fib in np.arange(a):
            diff = np.diff(wave[fib,:])
            diff_array = np.hstack([diff,diff[-1]])
            S[fib,:] = np.interp(std_wave, wave[fib,:], twi[fib,:]/ diff_array, 
                          left=-999., right=-999.)
        avg_spec_list.append(S)
    average = np.vstack(avg_spec_list)
    average = np.ma.array(average, mask=(average==-999.))
    average = biweight_location(average,axis=(0,))
    fiber_to_fiber = []
    for avg_spec, wave in zip(avg_spec_list, wave_list):
        y = func(avg_spec / average, (1,smooth_scale))
        fiber_to_fiber.append(y)
    return fiber_to_fiber, std_wave


def starextract(argsv=None, seeing=1.5, wavescale=1.9):
    '''
    Parameters
    ----------
    args : list of strings, optional
        command line
    '''
    log = setup_logging()
    parser = ap.ArgumentParser(description="""Calculate the fiber
                                     to fiber across a shot.""",
                                     parents=[astro_parent, ])
                                     
    parser.add_argument('--pos', nargs='?', type=str, 
                            help='''"x,y"''', default=None)
                            
                          
    args = parser.parse_args(argsv)
    if args.pos is None:
        log.error('No position argument "--pos" given.')
        sys.exit(1)
    args.ifuslots = args.ifuslots.replace(" ", "").split(',') 
    
    try:
        x,y = [float(i) for i in args.pos.replace(" ", "").split(',')] 
    except:
        log.error('Could not proper read in %s as an x,y position' %args.pos)
        sys.exit(1)
        
    for ifuslot in args.ifuslots:
        pass
        
def get_astrometry(argsv=None, collapse_lim1=5050, collapse_lim2=5400, 
                   imscale=0.5, seeing=1.5):
    '''
    Parameters
    ----------
    args : list of strings, optional
        command line
    '''
    log = setup_logging()
    parser = ap.ArgumentParser(description="""Calculate the fiber
                                     to fiber across a shot.""",
                                     parents=[astro_parent, ])
                                     
    parser.add_argument('--write_to_file', action="count", default=0,
                        help=("Add argument if you want to write images "
                        "to file"))       
                          
    args = parser.parse_args(argsv)
    args.ifuslots = args.ifuslots.replace(" ", "").split(',')
    
    # Get Shots 
    shot_list = get_shot_info_from_input(args)
    
    # Grab astrometry from headers and build projection from 1st shot
    astrom_list = get_astrometry_info_from_input(args, args.ifuslots[0])

    Ast = Astrometry(astrom_list[0][0]*15., astrom_list[0][1], 
                     astrom_list[0][2], fplane_file=args.fplane)
    
    # Calculate center of ifuslots and size and then grab images and catalogs
    # for that area from skyview/vizier                 
    ra, dec, size = get_image_limits(args, Ast)
    log.info('Getting archive image for RA: %0.5f, '
             'Dec: %0.5f, Size: %0.1f"'%(ra,dec,size))
    archive_image, archive_wcs = retrieve_image_SkyView(ra, dec, size, imscale)
    catalog = querySources(ra, dec, size)    
    cat = SkyCoord(catalog[:,0],catalog[:,1], frame='fk5', 
                   unit=(unit.degree, unit.degree))  
                   
    image_list, sources_list, matches_list, back_list = [],[],[],[]
    for shot, astrom in zip(shot_list, astrom_list):
        image_shot, sources_shot, matches_shot, back_shot = [],[],[],[]
        amp_list = []
        for k, ifuslot in enumerate(args.ifuslots):
            avg_scatter, std_scatter, amps = investigate_background(argsv, shot, 
                                                                    ifuslot)
            sel1 = np.abs(avg_scatter / astrom[4] * 360.) < 1.
            for amp in amps[~sel1]:
                log.info('Amp, %s, for IFUSLOT, %s, was rejected for high'
                         ' background.' %(amp,ifuslot))
            sel2 = std_scatter < 4.
            for amp in amps[~sel2]:
                log.info('Amp, %s, for IFUSLOT, %s, was rejected for significant'
                         ' scatter in the background.' %(amp,ifuslot))
            sel = np.where(sel1 * sel2)[0]
            amp_list.append(amps[sel])
        # Get Fiber to Fiber for shots to build collapsed images for each ifu
        FtF, std_wave = fiber_to_fiber_shot(shot, args.folder, args.instr, 
                                             args.ifuslots, amp_list)
        xl = np.searchsorted(std_wave, collapse_lim1, side='left')
        xh = np.searchsorted(std_wave, collapse_lim2, side='right') 
                                        
        nifus = len(args.ifuslots)
        
        for k, ifuslot in enumerate(args.ifuslots):
            # Making collapsed image for ifuslot
            log.info('Making collapsed image for %s on %s, %s, %s' 
                     %(ifuslot, shot[0], shot[1], shot[2]))

            if len(amp_list[k])==0:
                continue
            image = make_collapsed_image(args, shot, ifuslot, std_wave,
                                         FtF[k], imscale, seeing, xl, xh,
                                         amps=amp_list[k])            
                                       
            # Measure Background
            log.info('Measuring Background for %s on %s, %s, %s' 
                     %(ifuslot, shot[0], shot[1], shot[2]))
            sub, bkg = measure_image_background(image)
            
            # Detecting Sources in Sky-subtacted collapsed frame
            log.info('Detecting Sources in %s on %s, %s, %s' 
                     %(ifuslot, shot[0], shot[1], shot[2]))
            segm, sources = detect_in_image(sub, bkg, fwhm=seeing, 
                                            scale=imscale)
            
            # Make cutout of larger archive image/catalog for ifuslot
            
            
            log.info('Making plot for %s on %s, %s, %s' 
                     %(ifuslot, shot[0], shot[1], shot[2]))
            

            
            image_shot.append(sub)
            back_shot.append(bkg.background)
            sources_shot.append([sources, ifuslot])

        image_list.append(image_shot)
        sources_list.append(sources_shot)
        matches_list.append(matches_shot)
        back_list.append(back_shot)

        for i in np.arange(2):
            matches_shot = match_to_database(sources_shot, image_shot, cat, Ast, 
                                             imscale,thresh=8./(i+1))    
            mat = np.array(matches_shot)
            x,q = kmeans2(mat[:,6:8], 2-i)
            C = Counter(q)
            sel = np.where(q == C.most_common(1)[0][0])[0]
            dra = mat[sel,6]
            ddec = mat[sel,7]
            ast_keep = np.array([biweight_location(dra)*np.cos(np.deg2rad(astrom[1])),
                                 biweight_location(ddec),
                                 biweight_midvariance(dra)*np.cos(np.deg2rad(astrom[1])),
                                 biweight_midvariance(ddec)])
            print(ast_keep*3600.)
            Ast.dra = Ast.dra + biweight_location(dra)
            Ast.ddec = Ast.ddec + biweight_location(ddec)
            Ast.update_projection()
            plt.figure(figsize=(6,6))
            from scipy.stats import rankdata
            plt.scatter(mat[:,6]*3600.,mat[:,7]*3600., c=rankdata(mat[:,0],'dense'))
            plt.scatter(x[:,0]*3600.,x[:,1]*3600.,marker='x',color='k')
            plt.show()
        make_shot_plot(shot, nifus, args, Ast, cat, archive_image, archive_wcs, 
                       imscale, image_shot, sources_shot)
    return image_list, sources_list, matches_list , back_list, FtF


def make_shot_plot(shot, nifus, args, Ast, cat, archive_image,
                   archive_wcs, imscale, sub_list, 
                   source_list):
    fig = plt.figure(figsize=(4*2, 4*nifus))
    grid = [nifus, 2]
    for k, ifuslot in enumerate(args.ifuslots):

        sub = sub_list[k]
        sources = source_list[k]
        # Get projection for collapsed frame            
        Ast.get_ifuslot_projection(ifuslot, imscale, sub.shape[1]/2.,
                                       sub.shape[0]/2.)
        ra, dec = Ast.get_ifuslot_ra_dec(ifuslot)
        position = SkyCoord(ra, dec, unit="deg", frame='fk5')
        d2d = position.separation(cat)
        sel = np.where(d2d.arcsec < 40.)[0]
        size = 60./imscale
        cutout = Cutout2D(archive_image, position, (size,size), 
                          wcs=archive_wcs)
        make_image_subplot(fig, sub, Ast.tp.wcs, grid, 1, k,
                           use_projection=True, title=ifuslot, 
                           vval=[-1,5], use_norm=True, cmap='Greys')
        x,y = [],[]
        for source in sources[0]:
            x.append(source['xcentroid'])
            y.append(source['ycentroid'])
        if len(x):
            x = np.array(x)
            y = np.array(y)
            plt.scatter(x, y, color='none', edgecolor='r', marker='o', s=100)
        make_image_subplot(fig, cutout.data, cutout.wcs, [nifus, 2], 2, k,
                               vval=[-.05,.3],use_norm=True, cmap='Greys')  
        if len(x):
            x1, y1 = Ast.convert_ifuslot_xy_to_new_xy(x, y, cutout.wcs)
            plt.scatter((x1-cutout.data.shape[1]/2.+.5)
                        *cutout.wcs.wcs.cdelt[0]*3600., 
                        (y1-cutout.data.shape[0]/2.+.5)
                        *cutout.wcs.wcs.cdelt[1]*3600., color='none', 
                        edgecolor='r', marker='o', s=100)
        if len(sel):
            x2, y2 = cat[sel].to_pixel(cutout.wcs)
            plt.scatter((x2-cutout.data.shape[1]/2.+.5)
                        *cutout.wcs.wcs.cdelt[0]*3600., 
                        (y2-cutout.data.shape[0]/2.+.5)
                        *cutout.wcs.wcs.cdelt[1]*3600., color='g', 
                        marker='x', s=75)
    path = build_path(args.folder, args.instr, shot[0], shot[1], shot[2])
    fig.savefig(op.join(path, 'Collapsed_images_astrometry.pdf'), dpi=200)  
    plt.close(fig)    

def inspect_significance_map(image, trace, Sig, oversample_value, sigthresh, 
                             amps):

    # Get the arrays from the lists
        # Adjust trace for collected amps
    nfibers = trace.shape[0] / len(amps)
    D = image.shape[0] / len(amps)
    for i in np.arange(len(amps)):
        trace[(i*nfibers):((i+1)*nfibers),:] += D*i
            
    import pyds9
    ds9 = pyds9.DS9()
    ds9.set_np2arr(image)
    ds9.set('scale mode zscale')
    fibind, xind = np.where(Sig>sigthresh)
    fib_dict = {}
    for i,fib in enumerate(Sig):
        fib_dict[i] = 0
    for fib, xind in zip(fibind,xind):
        fib_dict[fib] +=1
        if fib_dict[fib] < 50:
            x = 1.*xind/oversample_value + 1
            y = trace[fib, int(x)-1]
            s = Sig[fib,xind]
            ds9.set('region','image; circle %f %f %f # color=red' %(x,y,s))
    
def calculate_significance(spectrum, wavelength, ftf, trace, error, image=None,
                           interactive=False, oversample_value=3, sigma=1.3, 
                           wave_step=4, sigthresh=4, cosmic_avoidance=4):
    '''
    Calculates the significance of map from extracted spectra.  
    '''
    # Oversample and Convolve
    convolved_spec = []
    oversampled_wave = []
    oversampled_ftf = []
    oversampled_spec = []
    for spec, wave, fiber_to_fiber in zip(spectrum, wavelength, ftf):
        s = np.diff(wave) / oversample_value
        t = np.hstack([s, s[-1]])
        wave_list = []
        for i in np.arange(oversample_value):
            wave_list.append([wave + i*t])
        oversampled_wave.append(np.sort(np.hstack(wave_list)).ravel())
        oversampled_spec.append(np.interp(oversampled_wave[-1], wave, 
                                     spec).ravel())
        oversampled_ftf.append(np.interp(oversampled_wave[-1], wave, 
                                     fiber_to_fiber).ravel())
        kernel = Gaussian1DKernel(stddev=sigma*oversample_value)
        convolved_spec.append(convolve(oversampled_spec[-1], kernel))
    spectrum = np.array(convolved_spec)
    wavelength = np.array(oversampled_wave)
    ftf = np.array(oversampled_ftf)
    
    # Calculate Noise
    wvmin = wavelength.min()
    wvmax = wavelength.max()
    wave_array = np.arange(wvmin,wvmax+2*wave_step, wave_step)
    wave_array_fit = wave_array[:-1]+wave_step/2.
    xind=[]
    yind=[]
    lx = 0
    ly = len(wave_array)-1
    for i in np.arange(ly):
        x,y = np.where((wavelength>=wave_array[i]) * (wavelength<wave_array[i+1]))
        lx = np.max([lx,len(x)])
        xind.append(x)
        yind.append(y)
    A = np.ones((ly,lx))*-999.
    for i in np.arange(ly):
        A[i,:len(xind[i])] = np.where(ftf[xind[i],yind[i]]>1e-8, 
                                      spectrum[xind[i],yind[i]]/ftf[xind[i],yind[i]],
                                      0.0)
    B = np.ma.array(A, mask=(A==-999.).astype(np.int))    
    C = biweight_midvariance(B,axis=(1,))
    V = np.zeros(spectrum.shape)
    for i in np.arange(ly):
        V[xind[i],yind[i]] = np.interp(wavelength[xind[i],yind[i]],wave_array_fit,C)
    
    # Calculate Significance
    Sig = spectrum/(V*ftf)
    Sig[~np.isfinite(Sig)] = -999.
    
    # Flag values near cosmics
    cyind, cxind = np.where(error==-1)
    
    for xind,yind in zip(cxind,cyind):
        trace_a = trace[:,xind]
        fibs = np.where(np.abs(trace_a-yind)<cosmic_avoidance)[0]
        for fib in fibs:
            lx = (xind-cosmic_avoidance)*oversample_value
            hx = (xind+cosmic_avoidance)*oversample_value + 1
            Sig[fib,lx:hx] = -999. 


    return Sig, wavelength
    
def make_collapsed_image(args, shot, ifuslot, std_wave, FtF, imscale,
                         seeing, xl, xh, amps=None, immin=-25., immax=25.):
    wave, spec, ifupos = read_in_multifiles(args.folder, shot[0], shot[1], 
                                   shot[2], ifuslot, args.instr, amps=amps,
                                   extensions=['wavelength',
                                               'spectrum',
                                               'ifupos'])
    # Get the arrays from the lists
    wave, spec, ifupos = wave[0], spec[0], ifupos[0]
    
    # Interpolate spectra to standard wavelength
    std_spec = np.zeros((wave.shape[0],len(std_wave)))   
    for fib, wv in enumerate(wave):
        std_spec[fib,:] = np.interp(std_wave, wv, spec[fib,:], 
                  left=-999., right=-999.)
                  
    # Correct for Fiber to Fiber variations
    corr_spec = np.where(FtF>1e-8, std_spec/FtF, -999.)
    corr_spec = np.ma.array(corr_spec, mask=(corr_spec==-999.))
    
    # Collapse Fibers
    fiber_vals = biweight_location(corr_spec[:,xl:xh],axis=(1,))
    
    # Create collapse fiber image
    x = np.arange(immin, 
                  immax+imscale, imscale)
    y = np.arange(immin, 
                  immax+imscale, imscale)
    xgrid, ygrid = np.meshgrid(x, y)
    d = np.zeros((ifupos.shape[0],)+xgrid.shape)
    w = np.zeros((ifupos.shape[0],)+xgrid.shape)
    farea = imscale**2 / (0.75**2*np.pi) 
    for i in xrange(len(x)):
        for j in xrange(len(y)):
            d[:,j,i]= np.sqrt((ifupos[:,0] - xgrid[j,i])**2 + 
                              (ifupos[:,1] - ygrid[j,i])**2)
            w[:,j,i] = np.exp(-1./2.*(d[:,j,i]/seeing)**2)
    ws = w.sum(axis=0)
    return ((fiber_vals[:,np.newaxis,np.newaxis]*w).sum(axis=0) / ws * farea)
  
def measure_image_background(image):
    '''
    image
    '''
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = BiweightLocationBackground()
    bkg = Background2D(image, (20, 20), filter_size=(3, 3),
                        sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    return image-bkg.background, bkg


def detect_in_image(image_sub, bkg, thresh=3., fwhm=1.5, scale=0.5):
    '''
    Make detections in sky-subtracted image
    image_sub, bkg
    '''
    threshold = (thresh * bkg.background_rms)
    fwhm_i = fwhm /scale
    sigma_i = fwhm_i * gaussian_fwhm_to_sigma
    kernel_size = int(sigma_i*4.)
    kernel = Gaussian2DKernel(sigma_i, x_size=kernel_size, y_size=kernel_size)
    kernel.normalize()
    segm = detect_sources(np.array(image_sub), threshold, npixels=5, 
                          filter_kernel=kernel)
    if segm.array.sum():
        segm = deblend_sources(image_sub, segm, npixels=5,
                                       filter_kernel=kernel)
    d = DAOStarFinder(fwhm=fwhm_i, threshold=0.15*thresh*bkg.background_rms_median,
                      exclude_border=True)
    sources = d(image_sub)
    return segm, sources

def make_image_subplot(fig, image, wcs, grid, i, k, title=None,
                       vval=None, use_norm=False, cmap='Greys',
                       use_projection=False):
    figkwargs={}
    kwargs = {}
    if use_projection:
        figkwargs['projection']=wcs
    else:
        extent = [-image.shape[1]/2.*wcs.wcs.cdelt[0]*3600.,
                  image.shape[1]/2.*wcs.wcs.cdelt[0]*3600.,
                  -image.shape[0]/2.*wcs.wcs.cdelt[1]*3600.,
                  image.shape[0]/2.*wcs.wcs.cdelt[1]*3600.]
        kwargs['extent'] = extent
    fig.add_subplot(grid[0], grid[1], i + (k*grid[1]), **figkwargs)
    kwargs['origin'] = 'lower'
    kwargs['interpolation'] = 'none'
    kwargs['cmap'] = cmap
    if use_norm:
        kwargs['norm'] = ImageNormalize(stretch=AsinhStretch())
    if vval is not None:
        kwargs['vmin'] = vval[0]
        kwargs['vmax'] = vval[1]
    plt.imshow(image, **kwargs)
    if title is not None:
        plt.title(title)
    ax = plt.gca()
    plt.xlim(ax.get_xlim())
    plt.ylim(ax.get_ylim())
    
def build_source_apertures(props, wcsi, wcse):
    r = 1.5    # approximate isophotal extent
    apertures = []
    for prop in props:
        ra, dec = wcsi.wcs_pix2world(prop.xcentroid.value, 
                                     prop.ycentroid.value,1)
        position = wcse.wcs_world2pix(ra,dec,1)
        a = prop.semimajor_axis_sigma.value * r
        b = prop.semiminor_axis_sigma.value * r
        theta = prop.orientation.value
        apertures.append(EllipticalAperture(position, a, b, theta=theta))
    return apertures

def call_skyview(survey, ra, dec, size, size_pix):
    try:
        paths = SkyView.get_images(position='%3.5f %2.5f' % (ra, dec),
                                   survey=survey, coordinates='J2000',
                                   height=size*aunit.arcsec, width=size*aunit.arcsec,
                                   pixels=str(size_pix))
    except:
        paths=[]
    return paths

def retrieve_image_SkyView(ra, dec, size, scale):
    log=setup_logging()
    size_pix = int(1.*size/scale)
    log.info('Retrieving image from SDSS')
    paths = call_skyview('SDSSr', ra, dec, size, size_pix)
                 
    if paths:
        wcs = WCS(paths[0][0])
        imarray = paths[0][0].data
    else:
        # Default to DSS2 Blue
        log.info('Retrieving image from SDSS did not work')
        log.info('Retrieving image from DSS2 Blue')
        paths = call_skyview('DSS2 Blue', ra, dec, size, size_pix)
        wcs = WCS(paths[0][0])
        imarray = paths[0][0].data

    return imarray, wcs

def get_image_limits(args, Ast):
    ra_dec = []
    for ifuslot in args.ifuslots:
        ra_dec.append(Ast.get_ifuslot_ra_dec(ifuslot))
    ra_dec = np.vstack(ra_dec)

    RA = np.mean(ra_dec[:,0])
    DEC = np.mean(ra_dec[:,1])
    tp = TP(RA, DEC, 0., 1., 1.)
    x,y = tp.raDec2xy(ra_dec[:,0],ra_dec[:,1])
    size = np.max([np.abs(x).max(), np.abs(y).max()])*2.+80.
    return RA, DEC, size

def match_to_database(source_list, image_list, cat, Ast, imscale, thresh=8.):
    '''
    RA, Decdetection catalogs and vizier (astroquery) match?  vo conesearch?    
    '''
    matches = []
    for sources, image in zip(source_list, image_list):
        ifuslot = sources[1]
        Ast.get_ifuslot_projection(ifuslot, imscale, image.shape[1]/2.,
                                       image.shape[0]/2.)
        for source in sources[0]:
            ra, dec = Ast.tp_ifuslot.wcs.wcs_pix2world(source['xcentroid'], 
                                                       source['ycentroid'],1)
            c = SkyCoord(ra, dec, unit=(unit.degree, unit.degree))
            idx, d2d, d3d = match_coordinates_sky(c, cat, nthneighbor=1) 
            if d2d.arcsec[0] < thresh:
                dra =  cat.ra[idx].deg - float(ra)
                ddec = cat.dec[idx].deg - float(dec)
                matches.append([int(ifuslot),int(idx), 
                                float(ra), float(dec), 
                                cat.ra[idx].deg, 
                                cat.dec[idx].deg, dra, ddec,
                                d2d.arcsec[0]])
    return matches
    
    
def querySources(ra, dec, boxsize, maxsources=10000):
    try:
        sources = querySDSS(ra, dec, boxsize, maxsources)
    except:
        sources = queryGAIA(ra, dec, boxsize, maxsources)
    if len(sources)==0:
        sources = queryGAIA(ra, dec, boxsize, maxsources)
    return sources

def querySDSS(ra, dec, boxsize, maxsources=10000): 
    """
    Queries USNO_A2.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. B R 0. 0.
    """
    log = setup_logging()
    log.info("querySDSScat: ra      = %f ", ra)
    log.info("querySDSScat: dec     = %f ", dec)
    log.info("querySDSScat: boxsize = %f ", boxsize)
    
    vquery = Vizier(columns=['gmag', 'RAJ2000', 'DEJ2000'], 
                    row_limit = maxsources) 
 
    field = coord.SkyCoord(ra=ra, dec=dec, 
                           unit=(unit.deg, unit.deg), 
                           frame='fk5')
    Data =  vquery.query_region(field, 
                               width=("%fd" % (boxsize/3600.)), 
                               catalog="II/294/sdss7")[0] 
    oo = []
    for i, obj in enumerate(Data['gmag']):
        ra = Data['RAJ2000'][i]
        dec = Data['DEJ2000'][i]
        B = Data['gmag'][i]
        if np.any([j for j in Data.mask[i]]):
            continue
        oo.append([ra, dec, B])
    return np.array(oo)

def queryGAIA(ra, dec, boxsize, maxsources=10000): 
    """
    Queries USNO_A2.
    ra  = center RA of field
    dec = center DEC of field
    radius = determines size of radius around ra/dec for which sources should
              be retrieved
    return array of stars with format
    IDa IDb RA DEC 0. B R 0. 0.
    """
    log = setup_logging()
    log.info("queryGAIAcat: ra      = %f ", ra)
    log.info("queryGAIAcat: dec     = %f ", dec)
    log.info("queryGAIAcat: boxsize = %f ", boxsize)
    
    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 
                             'phot_g_mean_mag','_RAJ2000', '_DEJ2000'], 
                    row_limit = maxsources) 
 
    field = coord.SkyCoord(ra=ra, dec=dec, 
                           unit=(unit.deg, unit.deg), 
                           frame='fk5')
    Data =  vquery.query_region(field, 
                               width=("%0.5fd" % (boxsize/3600.)), 
                               catalog="I/337/gaia")[0] 
    oo = []
    for i, obj in enumerate(Data['Source']):
        ra = Data['_RAJ2000'][i]
        dec = Data['_DEJ2000'][i]
        B = Data['__Gmag_'][i]
        if np.any([j for j in Data.mask[i]]):
            continue
        oo.append([ra, dec, B])
    return np.array(oo)

def investigate_background(args, shot, ifuslot):
    log = setup_logging()
    parser = ap.ArgumentParser(description="""Calculate the fiber
                                     to fiber across a shot.""",
                                     parents=[astro_parent, ])
    args = parser.parse_args(args)
    args.ifuslots = args.ifuslots.replace(" ", "").split(',')
    path = build_path(args.folder, args.instr, shot[0], shot[1], shot[2])
    fn_search = op.join(path,'multi_*_%s_*_*.fits' %(ifuslot))
    fn = glob.glob(fn_search)
    if not fn:
        log.error('No files found for %s' %fn_search)
        sys.exit(1)
    avg,std,amp = [],[],[]
    for f in fn:
        F = fits.open(f)
        image = F['image'].data
        trace = F['trace'].data
        N,D = image.shape
        Fibers = []
        for i,tr in enumerate(trace):
            Fibers.append(Fiber(D, i+1, F[0].header['RAWPATH'], 
                                F[0].header['RAWFN']))
            Fibers[-1].trace = tr
        get_indices(image, Fibers, F[0].header['FSIZE'])
        t = []
        ygrid,xgrid = np.indices(image.shape)
        ygrid = 1. * ygrid.ravel() / N
        xgrid = 1. * xgrid.ravel() / D
        image = image.ravel()
        s = np.arange(N*D)
        for fiber in Fibers:
            t.append(fiber.D*fiber.yind + fiber.xind)
        t = np.hstack(t)
        t = np.array(t, dtype=int)
        ind = np.setdiff1d(s,t)
        mask = np.zeros((N*D))
        mask[ind] = 1.
        sel = np.where(mask==1.)[0]
        avg.append(biweight_location(image[sel]))
        std.append(biweight_midvariance(image[sel]))
        amp.append(F[0].header['AMP'])
    return np.array(avg), np.array(std) , np.array(amp)

def throughput(): 
    pass

