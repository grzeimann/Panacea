# -*- coding: utf-8 -*-
"""
Fiber Utilities
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

from utils import biweight_location, biweight_midvariance, is_outlier
from utils import biweight_filter
from scipy.optimize import nnls
import scipy
from scipy.signal import medfilt, savgol_filter

from scipy.interpolate import interp1d, NearestNDInterpolator, LinearNDInterpolator
from scipy.linalg import lstsq
from numpy.polynomial.polynomial import polyvander2d
import matplotlib.pyplot as plt
import matplotlib
import os.path as op
import numpy as np
import time
import sys
from astropy.convolution import convolve, Gaussian1DKernel
from bspline import Bspline
import splinelab

def str2bool(v):
    '''
    Convert a string to a boolean.  If the string is in the list below,
    then True is returned.    
    
    '''
    
    return v.lower() in ("y", "yes", "true", "t", "1")

def get_float(answer, question, previous):
    '''
    Examines string return for various cases.
    If 'q' then quit.  If '' then use previous value.  If float, return it.
    
    '''
    try:
        value = float(answer)
    except ValueError:
        if answer =='q':
            sys.exit(1)
        if answer =='':
            value = previous
        else:
            answer = raw_input(question)
            get_float(answer, question, previous)
    return value

def get_int(answer, question, previous):
    '''
    Examines string return for various cases.
    If 'q' then quit.  If '' then use previous value.  If int, return it.
    
    '''
    try:
        value = int(answer)
    except ValueError:
        if answer =='q':
            sys.exit(1)
        elif answer =='':
            value = previous
        else:
            answer = raw_input(question)
            get_float(answer, question, previous)
    return value

def fit_scale_wave0(init_scale, init_wave0, xi, xe, D, sun_wave, ysun, data, 
                    fixscale=False, dofit=True, buff=10,
                    plot=False, res=1.9): 
    '''
    Fitting function for linear wavelength solution of a small column section
    of one fiber.  
    '''
    if fixscale:
        def f(params, sel1):
            wv = init_scale * np.arange(D) + params[0]
            model = np.interp(wv[xi:xe+1], sun_wave, ysun, left=0.0, right=0.0)
            return model[sel1] - data[sel1]
        if dofit:
            params0 = np.array([init_wave0])
            sel = is_outlier(f(params0,np.ones(data.shape,dtype=bool)))<1
            sol = scipy.optimize.leastsq(f, params0, args=(sel))[0]
            chi2 = f(sol, sel+True )**2
        else: 
            params0 = np.arange(init_wave0-buff, init_wave0+buff, res/10.)
            sel = np.ones(data.shape,dtype=bool)
            chi2_manual = np.zeros(params0.shape)
            for i,p in enumerate(params0):
                chi2_manual[i] = (f([p], sel)**2).sum() 
            if plot:
                plt.figure()
                plt.plot(params0,chi2_manual)
                plt.yscale('log')
                plt.show()
                raw_input("Please press enter")
                plt.close()
            sol = [params0[chi2_manual.argmin()]]
            chi2 = f(sol, sel)**2
    else:
        def f(params, sel1):
            wv = params[0] * np.arange(D) + params[1]
            model = np.interp(wv[xi:xe+1], sun_wave, ysun, left=0.0, right=0.0)
            return model[sel1] - data[sel1]
                
        params0 = np.array([init_scale, init_wave0])
        sel = is_outlier(f(params0, np.ones(data.shape,dtype=bool)))<1
        if np.sum(sel)>3:
            sol = scipy.optimize.leastsq(f, params0, args=(sel))[0]
            chi2 = f(sol, sel+True )**2
        else:
            sol = params0
            chi2 = 1e6*np.ones(data.shape)
    return sol, chi2

def find_maxima(x, y, y_window=3, interp_window=2.5, repeat_length=2,
                first_order_iter=5, max_to_min_dist=5):
    '''
    This function finds maxima related to fiber peaks.  The centroid
    is defined iteratively over an interpolated field.
    
    :param x:
        The row values for a column slice.
    :param y:
        The counts for a given column slice.
    :param y_window:
        In [x - y_window : x + y_window] evaluate the maximum value.
    :param interp_window:
        In [x - interp_window : x + interp_window] interpolate y to new grid
        and find the first order moment (centroid)
    :param repeat_length:
        In the moving window the maximum value has to repeat at least this 
        many times + 1 to keep it as a potential fiber centroid
    :param first_order_iter:
        The number of iterations to re-calculate the first order moment
        over the interpreted grid.
    :param max_to_min_dist:
        Initial maxima are found over [min, next min] grid and the min/next min
        is fixed to be at most max_to_min_dist away
    '''
    a = len(y)
    mxval = np.zeros((a-2-y_window*2+1,))
    mnval = np.zeros((a-2-y_window*2+1,))
    k=0
    for j in xrange(y_window,a - 1 - y_window):
        mxval[k] = np.argmax(y[(j-y_window):(j+y_window+1)])+j-y_window
        mnval[k] = np.argmin(y[(j-y_window):(j+y_window+1)])+j-y_window
        k+=1
    # Get unique values
    umnval = np.unique(mnval)
    umxval = np.unique(mxval)
    mnkeep = []
    mxkeep = []
    # Only keep those which are max/min for at least repeat_length
    for j in umnval:
        if (j==mnval).sum() > repeat_length:
           mnkeep.append(j)
    for j in umxval:
        if (j==mxval).sum() > repeat_length:
           mxkeep.append(j) 
           
    # Make lists arrays
    mnkeep = np.array(mnkeep,dtype=int)
    mxkeep = np.array(mxkeep,dtype=int)
    
    #refined peak values and heights
    peaks_refined = np.zeros(mxkeep.shape)
    peaks_height = np.zeros(mxkeep.shape)
    for j, val in enumerate(mxkeep):
        mn1 = np.max([int(val - max_to_min_dist),0])
        mn2 = np.min([int(val + max_to_min_dist+1),len(y)])
        p0 = np.polyfit(x[mn1:mn2],y[mn1:mn2],2)
        lwa = -p0[1] / (2.*p0[0])
        if np.isnan(lwa):
            lwa = ((y[mn1:mn2]*x[mn1:mn2]).sum()/(y[mn1:mn2]).sum())
        if lwa<mn1 or lwa>mn2:
            lwa = ((y[mn1:mn2]*x[mn1:mn2]).sum()/(y[mn1:mn2]).sum())

        # Iterating and interpolating to refine the centroid
        for k in xrange(first_order_iter):
            xp = np.linspace(lwa-interp_window, lwa+interp_window, num=50)
            yp = np.interp(xp, x[mn1:mn2], y[mn1:mn2], left=0, right=0)
            try:
                p0 = np.polyfit(xp,yp,2)
                lwa = -p0[1] / (2.*p0[0])
                if np.isnan(lwa):
                    lwa = (yp*xp).sum()/yp.sum()
                if lwa<mn1 or lwa>mn2:
                    lwa = (yp*xp).sum()/yp.sum()
            except:
                lwa = (yp*xp).sum()/yp.sum()
        peaks_refined[j] = lwa
        peaks_height[j] = y[val]
    return peaks_refined, peaks_height 

def calculate_wavelength(x, y, solar_peaks, solar_spec, window_size=80.,
                         isolated_distance=15., min_isolated_distance=4., 
                         init_lims=None, smooth_length=51, order=3, 
                         init_sol=None, height_clip=1.02, debug=False,
                         interactive=False, constrained_to_initial=False,
                         maxwavediff=5.):
    '''
    Outdated wavelength calculation using peaks in the normalized 
    twi/solar spectrum.
    
    '''
    if debug:
        t1 = time.time()
    if init_lims is None:
        init_lims = [np.min(solar_peaks[:,0]), np.max(solar_peaks[:,0])]
    # smooth_length has to be odd
    if smooth_length%2 == 0:
        smooth_length+=1
    yp = biweight_filter(y, smooth_length) / y
    p_loc, p_height = find_maxima(x, yp)
    ind_sort = np.argsort(solar_peaks[:,1])[::-1]
    matches = []
    if init_sol is None:
        init_sol = np.array([1.*(init_lims[1]-init_lims[0])/len(x), 
                             init_lims[0]])
        change_thresh=1
    else:
        change_thresh=999
        if constrained_to_initial:
            init_wave = np.polyval(init_sol,x)
    for i, ind in enumerate(ind_sort[:]):
        if (np.min(np.abs(solar_peaks[ind,0]-np.delete(solar_peaks[:,0],ind)))
                          <isolated_distance):
            print('Skipped due to proximity: %0.1f' 
                        %np.min(np.abs(solar_peaks[ind,0]
                        -np.delete(solar_peaks[:,0],ind))))
            print(isolated_distance)
            continue
        wv = np.polyval(init_sol, x)
        wvi = solar_peaks[ind,0]
        p_loc_wv = np.interp(p_loc, x, wv)
        selp = np.where((p_loc_wv > (wvi-window_size)) 
                        * (p_loc_wv < (wvi+window_size))
                        * (p_height > height_clip))[0]
        if selp.size:
            if interactive:
                xl = np.searchsorted(solar_spec[:,0],wvi-window_size)
                xu = np.searchsorted(solar_spec[:,0],wvi+window_size)
                ax = plt.axes([0.1,0.1,0.8,0.8])
                ax.step(solar_spec[xl:xu,0],solar_spec[xl:xu,1],where='mid')
                ax.step(wv,yp,'r-',where='mid')
                ax.scatter(wvi, solar_peaks[ind,1])
                #psel = ((solar_peaks[:,0]>wvi-window_size)*
                #        (solar_peaks[:,0]<wvi+window_size))
                #ax.scatter(solar_peaks[psel,0],solar_peaks[psel,1])
                ax.scatter(p_loc_wv[selp],p_height[selp],c='r')
                for s in selp:
                    ax.text(p_loc_wv[s],p_height[s]+.03,"%0.2f" %p_loc[s])
                ax.set_xlim([wvi-window_size,wvi+window_size])
                mn = solar_spec[xl:xu,1].min()
                mx = np.max([solar_spec[xl:xu,1].max(),np.max(p_height[selp])])
                rn = mx - mn
                ax.set_ylim([-.1*rn + mn, 1.1*rn+mn])
                plt.show()
                answer = raw_input("What is the xpos for wavelength {:0.0f}?".format(wvi))
                plt.close()
                try:
                    xv = float(answer)
                    matches.append([xv, wvi])
                except ValueError:
                    if answer=='off':
                        interactive=False
                    elif answer=='q':
                        sys.exit(1)
                    else:
                        continue                    
            if not interactive:
                s_ind = np.argmin(np.abs(p_loc_wv[selp]-wvi))
                m_ind = selp[s_ind]
                if len(selp)>1:
                    if (np.min(np.abs(p_loc[m_ind]-np.delete(p_loc[selp],s_ind)))
                              <isolated_distance):
                        continue
                xv = p_loc[m_ind]
                matches.append([xv, wvi])
        matches_array = np.array(matches)
        if len(matches)>change_thresh:
            if ((np.max(matches_array[:,0])
                    -np.min(matches_array[:,0]))/len(x) > 0.7):
                order_fit = np.min([np.max([len(matches)-2, 1]),order])
            else:
                order_fit = 1
            test_sol = np.polyfit(matches_array[:,0], matches_array[:,1],
                                  order_fit)
            if constrained_to_initial:                      
                newwave = np.polyval(test_sol,x)
                if np.max(newwave-init_wave)<maxwavediff:
                    init_sol = test_sol
            else:
                init_sol = test_sol

    bad = is_outlier(np.polyval(init_sol,matches_array[:,0])
                        -matches_array[:,1])
    print(np.sum(bad<1),np.sum(bad==1))
    order_fit = np.min([np.max([np.sum(bad<1)-2, 1]),order])
    init_sol = np.polyfit(matches_array[bad<1,0], matches_array[bad<1,1],
                              order_fit)
    std = biweight_midvariance(np.polyval(init_sol,
                                                  matches_array[bad<1,0])
                                       -matches_array[bad<1,1])
    print("Standard Deviation Offset: %0.2f" %std)
    if debug:
        t2=time.time()
        print("Time to finish wavelength solution for a single fiber: %0.2f" %(t2-t1))
    return np.polyval(init_sol,x), init_sol


def calculate_wavelength_new(x, solar_spec, fibers, fibn, group,
                              smooth_length=21, init_lims=None, order=3, 
                              init_sol=None, debug=False, interactive=False, 
                              nbins=21, wavebuff=100, plotbuff=85, 
                              fixscale=True, use_leastsq=False, res=1.9):
    L = len(x)
    if init_lims is None:
        init_lims = [np.min(solar_spec[:,0]), np.max(solar_spec[:,0])]
    if init_sol is not None:
        init_wave_sol = np.polyval(init_sol, 1. * x / L)
    y_sun = solar_spec[:,1]  
    lowfib = np.max([0,fibn-group])
    highfib = np.min([len(fibers)-1,fibn+group])
    y = np.array([biweight_filter(fibers[i].spectrum, smooth_length)
                  / fibers[i].spectrum 
                  for i in xrange(lowfib,highfib)])
    y = biweight_location(y,axis=(0,))
    bins = np.linspace(init_lims[0], init_lims[1], nbins)
    bins = bins[1:-1]
    scale = 1.*(init_lims[1] - init_lims[0])/L
    wv0 = init_lims[0]
    wave0_save = [] 
    scale_save = []
    x_save = []
    wave_save = []


def calculate_wavelength_chi2(x, solar_spec, fibers, fibn, group,
                              smooth_length=21, init_lims=None, order=3, 
                              init_sol=None, debug=False, interactive=False, 
                              nbins=21, wavebuff=100, plotbuff=85, 
                              fixscale=True, use_leastsq=False, res=1.9):
    '''
    Use a linear solution of the wavelength via a chi^2 fit from normalized
    twighlight spectrum and that of the template.  
    
    :param x:
        Column values
    :param y:
        Counts
    :param solar_spec:
        Two column vector array with wavelength and normalized flux of the 
        template spectrum
    :param smooth_length:
        Biweight smoothing length of the twi spectrum used for normalization.
    :param init_lims:
        Two entry list or tuple of the wavelength guess for the first and last
        column.  This is used as an initial linear solution.
    :param order:
        Polynomial order for the final wavelength solution of the fiber
    :param init_sol:
        Initial solution (polynomial) to start with instead of a linear guess
        from init_lims.
    :param debug:
        Used for random debug purposes, especially timing
    :param interactive:
        If true, the linear fits are interactive with feedback for slope and
        intercept.  This is helpful for debugging and difficult fits
    :param nbins:
        Number of linear bins fit for the ultimate polynomial solution
    :param wavebuff:
        The linear bins don't start at the first column but instead start and
        end at init_lims[0] + wavebuff and init_lims[1] - wavebuff.  This is
        critical for LRS2 where the throughput is very low in the blue side
    :param plotbuff:
        In the interactive mode, the plotted bin has wavelength limits that
        are +/-plotbuff the fit limits
    :param fixscale:
        Only fit the intercept in the linear fit of each bin
    '''
    L = len(x)
    if init_lims is None:
        init_lims = [np.min(solar_spec[:,0]), np.max(solar_spec[:,0])]
    if init_sol is not None:
        init_wave_sol = np.polyval(init_sol, 1. * x / L)
    y_sun = solar_spec[:,1]  
    lowfib = np.max([0,fibn-group])
    highfib = np.min([len(fibers)-1,fibn+group])
    y = np.array([biweight_filter(fibers[i].spectrum, smooth_length)
                  / fibers[i].spectrum 
                  for i in xrange(lowfib,highfib)])
    y = biweight_location(y,axis=(0,))
    bins = np.linspace(init_lims[0], init_lims[1], nbins)
    bins = bins[1:-1]
    scale = 1.*(init_lims[1] - init_lims[0])/L
    wv0 = init_lims[0]
    wave0_save = [] 
    scale_save = []
    x_save = []
    wave_save = []
    for j in xrange(len(bins)-1):
        if init_sol is not None:
            xl = np.searchsorted(init_wave_sol, bins[j])
            xh = np.searchsorted(init_wave_sol, bins[j+1])
            p0 = np.polyfit(x[xl:xh+1], init_wave_sol[xl:xh+1], 1)
            scale = p0[0]
            wv0 = p0[1]
        wv = np.arange(L)*scale + wv0
        xi = np.searchsorted(wv,bins[j])
        xe = np.searchsorted(wv,bins[j+1])
        swlow = bins[j]-wavebuff  
        swup = bins[j+1]+wavebuff
        xsl = np.searchsorted(solar_spec[:,0],swlow)
        xsu = np.searchsorted(solar_spec[:,0],swup)
        sun_wave = solar_spec[xsl:xsu,0]
        ysun = y_sun[xsl:xsu]
        flag = np.ones(y[xi:xe+1].shape, dtype=bool)
        content = False
        save=True
        while content is False:
            if init_sol is None and j==0:
                if use_leastsq:
                    xwv0 = np.arange(wv0-25,wv0+25,5)
                    chi2r = np.zeros(xwv0.shape)
                    solk = []
                    for cnt,wv0g in enumerate(xwv0):
                        sol, chi2 = fit_scale_wave0(scale, wv0g, xi, xe, len(x), sun_wave, ysun, y[xi:xe+1],
                                      fixscale=fixscale)
                        solk.append(sol)
                        if fixscale:
                            wv = np.arange(L) * scale + sol[0]
                        else:
                            wv = np.arange(L)*sol[0] + sol[1]           
                        model = np.interp(wv[xi:xe+1], sun_wave, ysun, left=0.0, right=0.0)
                        flag = is_outlier(model - y[xi:xe+1])
                        chi2r[cnt] = np.nansum(chi2[flag<1]) / (np.sum(flag<1)+1-2)
                    sol = solk[np.argmin(chi2r)]
                else:
                    sol, chi2 = fit_scale_wave0(scale, wv0, xi, xe, len(x), sun_wave, ysun, y[xi:xe+1],
                                                fixscale=fixscale, dofit=False, 
                                                buff=res*20., plot=False, res=res)  
            else:
                if use_leastsq:
                    sol, chi2 = fit_scale_wave0(scale, wv0, xi, xe, len(x), sun_wave, ysun, y[xi:xe+1],
                                  fixscale=fixscale)
                else:
                    sol, chi2 = fit_scale_wave0(scale, wv0, xi, xe, len(x), sun_wave, ysun, y[xi:xe+1],
                                                fixscale=fixscale, dofit=False, 
                                                buff=4., res=res)
            if fixscale:
                wv = np.arange(L) * scale + sol[0] 
            else:
                wv = np.arange(L)*sol[0] + sol[1]                                   
            model = np.interp(wv[xi:xe+1], sun_wave, ysun, left=0.0, right=0.0)
            flag = is_outlier(model - y[xi:xe+1]) < 1
            std = biweight_midvariance(model-y[xi:xe+1])
            if (std*4)>model.max():
                save=False
            if save:
                if fixscale:
                    wv0 = sol[0]
                else:
                    scale = sol[0]
                    wv0 = sol[1] 
                    

            if interactive: 
                wv = np.arange(L)*scale + wv0
                plt.figure(figsize=(8,6))
                ax = plt.axes([0.1,0.1,0.8,0.8])
                ax.step(sun_wave, ysun, where='mid')
                ax.step(wv[xi:xe+1],y[xi:xe+1],'r-',where='mid')
                ax.scatter(wv[xi:xe+1][flag<1], y[xi:xe+1][flag<1],marker='x',color='b')
                yi = np.interp(wv[xi:xe+1], sun_wave, ysun, left=0.0, right=0.0)
                ax.plot(wv[xi:xe+1], yi,'k-')
                ax.set_xlim([swlow+plotbuff, swup-plotbuff])
                mn = np.max([ysun.min(),y[xi:xe+1].min()])
                mx = np.min([ysun.max(),y[xi:xe+1].max()])
                rn = mx - mn
                ax.set_ylim([-1.0*rn + mn, 2.0*rn+mn])
                ax.text(swlow+plotbuff+10,1.0*rn+mn, "STD: %0.5f" %std)
                plt.show()
                answer = raw_input("Are you happy with the fit [%0.3f, %0.2f]?" %(scale, wv0))
                if answer == 'q':
                    sys.exit(1)
                if answer == 'Y':
                    interactive=False
                content = str2bool(answer)
                if content:
                    plt.close()
                    continue
                question = ("What is the scale" 
                                "(Current is %0.3f)?" % scale)
                answer = raw_input(question)
                scale = get_float(answer, question, scale)            
                question = ("What wavelength zeropoint" 
                                "(Current is %0.2f)?" % wv0)
                answer = raw_input(question)
                wv0 = get_float(answer, question, wv0)
                plt.close()
            else:
                content=True
        if save:
            x_save.append(x[xi:xe+1])
            wave_save.append(wv[xi:xe+1])
            wave0_save.append(wv0)
            scale_save.append(scale)
    
    init_sol = np.polyfit(1.*np.hstack(x_save)/L,np.hstack(wave_save), order)
    return np.polyval(init_sol,1. * x / L), init_sol        

def get_trace_from_image(image, x_step=4, y_window=3, x_window=5, repeat_length=2,
                         order=3, mx_cut=0.1, max_to_min_dist=5., debug=False,
                         first_order_iter=5, interp_window=3., xc=None):
    '''
    :param image:
        Image [FITS] from which to fit the trace.
    :param y_window:
        Within [-y_window, +y_window], we evaluate the index of the maximum
        value.  If that index is repeated :repeat_length: times then it is 
        logged as a potential maximum, unless it is already logged in which
        case it is skipped as to not repeat.
    :param x_window:
        We average the columns over [-x_window,+x_window] to avoid pixel 
        defects and other related issues.
    :param repeat_length:
        As discussed in :param y_window:, repeat_length is used to identify
        maxima via being the maximum point within [-y_window,+y_window] for a
        repeated length of :repeat_length:.
    :param order:
        The order of the polynomial used to smooth the trace via fitting.
    :param mx_cut:
        The cut used to define maxima from noise.  Must be greater than 
        mx_cut * <avg height of all peaks>.
    :param max_to_min_dist:
        Initial maxima are found over [min, next min] grid and the min/next min
        is fixed to be at most max_to_min_dist away
    :param debug:
        Used for a variety of purposes including timing
    :param interp_window:
        In [x - interp_window : x + interp_window] interpolate y to new grid
        and find the first order moment (centroid)
    :param first_order_iter:
        The number of iterations to re-calculate the first order moment
        over the interpreted grid        
    '''
    allfibers=[]
    a, b = image.shape
    x = np.arange(a)
    if xc is None:
        xc = np.arange(b)[x_window:(b-1-x_window):x_step]
    if debug:
        t1 = time.time()
    for i in xc:
        # Average 2*x_windwow columns at a time
        y = biweight_location(image[:,(i-x_window):(i+x_window+1)], axis=(1,))
        
        # Set empty arrays for maxima and minima found in 2*y_window searches
        mxval = np.zeros((a,))
        mnval = np.zeros((a,))
        k=0
        for j in xrange(a):
            mind1 = int(np.max([0, j-y_window]))
            mind2 = int(np.min([a, j+y_window]))
            mxval[k] = np.argmax(y[mind1:mind2])+mind1
            mnval[k] = np.argmin(y[mind1:mind2])+mind1
            k+=1
        # Get unique values
        umnval = np.unique(mnval)
        umxval = np.unique(mxval)
        mnkeep = []
        mxkeep = []
        # Only keep those which are max/min for at least repeat_length
        for j in umnval:
            if (j==mnval).sum() > repeat_length:
               mnkeep.append(j)
        for j in umxval:
            if (j==mxval).sum() > repeat_length:
               mxkeep.append(j) 
               
        # Make lists arrays
        mnkeep = np.array(mnkeep,dtype=int)
        mxkeep = np.array(mxkeep,dtype=int)
        
        #refined peak values and heights
        peaks_refined = np.zeros(mxkeep.shape)
        peaks_height = np.zeros(mxkeep.shape)

        # Loop through maxima and find the light-weighted average iteratively
        # Each iteration defines the new centroid and we interpolate the data
        # to a finer grid for the light-weighted average in a 2*y_window region
        for j, val in enumerate(mxkeep):
            ind = np.searchsorted(mnkeep,val)
            if ind==len(mnkeep):
                ind-=1
                mn1 = int(mnkeep[ind])
                mn2 = int(len(y))
            elif ind==0:
                mn1 = int(0)
                mn2 = int(mnkeep[ind]+1)
            else:
                mn1 = int(mnkeep[ind-1])
                mn2 = int(mnkeep[ind]+1)
            if val - mn1 > max_to_min_dist:
                mn1 = int(val - max_to_min_dist)
            if mn2 - val > (max_to_min_dist + 1):
                mn2 = int(val + max_to_min_dist + 1)
            lwa = ((y[mn1:mn2]*x[mn1:mn2]).sum()/(y[mn1:mn2]).sum())
            # Iterating and interpolating to refine the centroid
            for k in xrange(first_order_iter):
                xp = np.linspace(lwa-interp_window, lwa+interp_window, num=50)
                yp = np.interp(xp, x[mn1:mn2], y[mn1:mn2],left=0,right=0)
                lwa = (yp*xp).sum()/yp.sum()
            peaks_refined[j] = lwa
            peaks_height[j] = y[val]
        mh = biweight_location(peaks_height)
        peaks_final = peaks_refined[peaks_height > (mx_cut * mh)]
        allfibers.append(peaks_final[np.isfinite(peaks_final)]) 
    if debug:
        t2 = time.time()
        print("Time Taken for Trace: %0.2f" %(t2-t1))
    return allfibers, xc


def build_powerlaw(ygrid, Fibers, plaw_coeff, cols=None, spectrum=False):
    
    P = np.zeros(ygrid.shape+(len(Fibers),))
    
    xp = np.linspace(ygrid.shape[0]*-1, ygrid.shape[0], 2*ygrid.shape[0]*5)
    plaw = plaw_coeff[0] / (plaw_coeff[1] + plaw_coeff[2]
              * np.power(abs(xp / 2.5), plaw_coeff[3]))
    I = interp1d(xp, plaw, bounds_error=False, fill_value=(0,0), kind='linear')
    if cols is None:
        cols = np.arange(ygrid.shape[1], dtype=int)
    for i,fiber in enumerate(Fibers):
        if spectrum:
            P[:,cols,i] = I(ygrid[:,cols]-fiber.trace[np.newaxis,cols])*fiber.spectrum[np.newaxis,cols]
        else:
            P[:,cols,i] = I(ygrid[:,cols]-fiber.trace[np.newaxis,cols])
    return P
    
def fit_fibermodel_nonparametric(image, Fibers, plot=False, fsize=8., 
                                 fiber_group=4, bins=15, col_group=48,
                                 use_default=False,
                                 outfolder=None, sigma=2.5, power=2.5,
                                 plaw_coeff1=0.0004):
    '''
    This function orchestrates the calls to the sub-function, 
    fit_fibermodel_nonparametric_bins, and grids the image to fit the fiber 
    model in an empirical method.  The fibermodel is fit in column blocks
    but done for each fiber with a window of "fiber_group" size.  In this 
    grid, a single fibermodel is fit and assigned to the fiber of interest and
    central column of the group.
    
    :param image:
        Amplifier image
    :param Fibers:
        List of fiber class for each fiber in the amplifier
    :param plot:
        If true, individual plots will be produced for each fit
    :param fsize:
        Region in which fiber model is defined and at the end set to zero.
        The region is [-fsize, fsize] in pixels.
    :param fiber_group:
        Defines the number of fibers fit simultaneously for the fibermodel.
        In practice, fiber_group+1 fibers are fit and only half of the first
        and the last fiber in the group are used for convergence.
    :param bins:
        The initial number of bins used to describe the 
        empirical/non-parametric profile.  Two other bins are added near
        the peak for common asymetries that are in the data but not in the
        initial_profile that defines the x-values of the bins.  This assures
        more accurate modeling of the peak of the profiles.
    :param col_group:
        Number of columns used to fit a single profile.  This is not a moving
        window but done once in a block and moved to the next block.
    :param debug:
        All-purpose boolean for investigating issues.
    :param use_default:
        This skips fitting a profile and uses just a default empirical 
        profile with amplitudes defined by a gauss-hermite exponential with 
        sigma and power as the only input parameters.
    :param outfolder:
        This is used for the location of the individual plots.
    :param sigma:
        Parameter for initial profile.
        profile = np.exp(-1./2.*((np.abs(x/sigma))**power))
    :param power:
        Parameter for initial profile.
        
    '''
    plaw_coeff = np.array([plaw_coeff1, 0.5, 0.15, 1.0])
    a,b = image.shape 
    ygrid,xgrid = np.indices(image.shape)                       
    nfibs = len(Fibers) 
    ncols = b / col_group
    so = np.zeros((nfibs, ncols, bins+2))
    xcols = np.arange(col_group/2, int((ncols-1/2.)*col_group)+1, col_group)
    bins1, binx, sol = init_fibermodel(fsize=fsize, bins=bins, sigma=sigma,
                                       power=power)
    PL = build_powerlaw(ygrid, Fibers, plaw_coeff)
            
    for j in xrange(nfibs):
        for i in xrange(ncols):
            if not use_default:
                sol = fit_fibermodel_nonparametric_bins(image, xgrid, ygrid, 
                                                    Fibers, PL, fib=j, 
                                                    group=fiber_group, 
                                                    bins=bins1,
                                                    xlow=i*col_group, 
                                                    xhigh=(i+1)*col_group, 
                                                    plot=plot,
                                                    outfolder=outfolder,
                                                    sol=sol, binx=binx)                
            so[j,i,:] = sol         
    return so, xcols, binx

    
def init_fibermodel(fsize, bins, sigma=2.5, power=2.5):
    '''
    Initial fiber model defined using a modified gaussian.
    The second derivative defines the number of bins used.
    The number of bins are buffered by 0.0 on either side 
    indicating the boundary of the profile and that the profile goes to zero
    at a distance of fsize.
    
    :param fsize:
        The fiber model is defined in the y-direction 
        over the interval [-fsize, fsize]
    :param sigma:
        sigma for the modified gaussian
    :param power:
        instead of a power of 2, this is a free parameter for the gaussian
    '''
    # Define the bins   
    xt = np.linspace(-1*fsize,fsize,501.)
    yt = np.exp(-1./2.*((np.abs(xt/sigma))**power))
    yt /= np.sum(yt*(xt[2]-xt[1]))
    cf = np.cumsum(np.abs(np.diff(np.diff(yt))))
    cf /= np.max(cf)
    xtp = xt[1:-1]+(xt[1]-xt[0])/2.
    binx = np.interp(np.linspace(0,1,bins),cf,xtp)
    binx = np.hstack([binx[:(bins/2)], binx[bins/2]/2.+binx[bins/2-1]/2.,
                      binx[bins/2], binx[bins/2]/2.+binx[bins/2+1]/2.,
                      binx[((bins/2)+1):]])
    bins+=2
    sol = np.interp(binx,xt,yt,left=0.0,right=0.0)
    sol[0] = 0.0
    sol[-1] = 0.0
    return bins, binx, sol

def gauss_power(x, sig, n, norm=1.):
    return norm/np.sqrt(2*np.pi*sig**2)*np.exp(-np.abs(x)**n/(2*sig**n)) 

def power_law(x, c1, c2=.5, c3=.15, c4=1., sig=2.5):
        return c1 / (c2 + c3 * np.power(abs(x / sig), c4))

def get_norm_for_gauss_power(siglow=1.,sighigh=4.,powlow=1.5,powhigh=4.):
    S,P = np.meshgrid(np.linspace(siglow,sighigh,31), 
                      np.linspace(powlow,powhigh,26))
    x = np.linspace(-20,20,41)
    
    S = S.ravel()
    P = P.ravel()
    v = np.zeros(S.shape)
    for i,s in enumerate(S):
        v[i] = gauss_power(x, s, P[i]).sum()
    X = np.zeros((len(S),2))
    X[:,0] = S
    X[:,1] = P
    return LinearNDInterpolator(X,v, fill_value=1.)

def get_norm_for_power_law(c1, sig, xlen):
    if xlen%2:
        xlen+=1
    x = np.linspace(-xlen/2,xlen/2,xlen+1)
    return power_law(x,c1,sig=sig).sum()
    
def subtract_plaw(image, fibers):
    y = 1.*image
    neigh_array = np.arange(len(fibers))
    for neigh in neigh_array:
        yind = fibers[neigh].yind
        xind = fibers[neigh].xind
        fibm = fibers[neigh].plaw
        y[yind,xind] = image[yind,xind] - fibm*fibers[neigh].spectrum[xind]
    return y
    
def subtract_neighbor_fibers(image, fibers, fib, neighbor): 
    y = 1.*image
    neigh_array = np.delete(np.arange(-neighbor,neighbor+1), neighbor)
    for neigh in neigh_array:
        if (fib+neigh>=0) and ((fib+neigh)<len(fibers)):
            yind = fibers[fib+neigh].yind
            xind = fibers[fib+neigh].xind
            fibm = fibers[fib+neigh].core 
            y[yind,xind] = y[yind,xind] - fibm*fibers[fib+neigh].spectrum[xind]
    fibers[fib].subtracted = y[fibers[fib].yind, fibers[fib].xind]

def get_interp_list(fsize, bins, interp_kind):
    xinterp = np.linspace(-fsize,fsize, bins)
    lx = len(xinterp)
    f = np.zeros((lx,))
    interp_list = []
    for i in xrange(lx):
        f[i] = 1.0
        interp_list.append(interp1d(xinterp, f, bounds_error=False, 
                                    fill_value=(0,0), kind=interp_kind))
        f[i] = 0.0
    return interp_list
    
def remeasure_core(image, fibers, fib, neighbor, col_lims, bins, cols,
                   fsize, c1, c2, break_fix, interp_list):
    
    # Smarter binning practices  
    a,b = image.shape                 
    neigh_array = np.arange(-neighbor,neighbor+1)
    x=[]
    kind=[]
    kind_full = []
    col_select = cols[(cols>=col_lims[0])*(cols<col_lims[1])]
    fmult = int(fsize*2)
    cnt = 0
    sel_list = []
    for neigh in neigh_array:
        if (fib+neigh>=0) and ((fib+neigh)<len(fibers)):
            if len(fibers[fib+neigh].xind)<(b*fmult):
                sel = np.where(np.in1d(fibers[fib+neigh].xind, col_select))[0]
            else:
                sel = ((col_select * fmult)[:,np.newaxis] 
                        + np.arange(fmult)).ravel()
            if neigh==neigh_array[0] or fib+neigh==0:
                sel1 = np.where((fibers[fib+neigh].yoff>=0))[0]
                sel2 = np.intersect1d(sel, sel1)
            elif neigh==neigh_array[-1] or fib+neigh==(len(fibers)-1):
                sel1 = np.where((fibers[fib+neigh].yoff<=0))[0]
                sel2 = np.intersect1d(sel, sel1)
            else:
                sel2 = sel
            x.append(fibers[fib+neigh].yoff[sel])
            yind = fibers[fib+neigh].yind[sel2]
            xind = fibers[fib+neigh].xind[sel2] 
            kind.append(yind * b + xind) 
            kind_full = np.union1d(kind_full, kind[-1])
            sel_list.append([sel,sel2])
            cnt+=1
    kind_full = kind_full.astype(int)
    xinterp = np.linspace(-fsize,fsize, bins)
    lx = len(xinterp)
    
    
    weight_image = np.zeros((lx,a*b))
    cnt=0
    for neigh in neigh_array:
        if (fib+neigh>=0) and ((fib+neigh)<len(fibers)):
            sel = sel_list[cnt][0]
            sel2 = sel_list[cnt][1]
            yind = fibers[fib+neigh].yind[sel2]
            xind = fibers[fib+neigh].xind[sel2] 
            for i, v in enumerate(interp_list):
                weight_image[i,kind[cnt]] += (v(fibers[fib+neigh].yoff[sel2])
                                             *fibers[fib+neigh].spectrum[xind])
            cnt+=1
    Y = weight_image[:,kind_full].swapaxes(0,1)
    I = image.ravel()[kind_full]
    
    virus_sel = (xinterp<-break_fix) + (xinterp>break_fix)
    virus_sel1 = (xinterp<-break_fix)
    virus_sel2 = (xinterp>break_fix)
    p1 = [c1, c2+c1*fsize]
    p2 = [-c1,c2+c1*fsize]
    soli = np.hstack([np.polyval(p1,xinterp[virus_sel1]),
                      np.polyval(p2,xinterp[virus_sel2])])
    sel_v = np.where(virus_sel)[0]
    sel_s = np.where(~virus_sel)[0]
    for i,v in enumerate(sel_v):        
        I -= Y[:,v] * soli[i]     
    try:
        sol = nnls(Y[:,~virus_sel], I)[0]  
    except:
        print(col_lims, fib)
        sys.exit(1)
        return None
    solb = np.zeros((lx,))
    solb[sel_v] = soli
    solb[sel_s] = sol
    sol = solb
    return xinterp, sol
    
def new_fast_norm(image, Fibers, cols=None, mask=None):
    if mask is None:
        mask = np.zeros(image.shape)
    a,b = image.shape
    if cols is None:
        cols = np.arange(b)
    init_model = np.zeros((a,len(Fibers)))
    norm = np.zeros((len(Fibers),b))
    for col in cols:
        for i,fiber in enumerate(Fibers):
            xsel = np.where(fiber.xind==col)[0]
            init_model[fiber.yind[xsel],i] = fiber.core[xsel]
        if (mask[:,col]==0).sum()>(a*1./2.):
            norm[:,col] = lstsq(init_model[mask[:,col]==0,:], 
                                image[mask[:,col]==0,col])[0]
                            
    for i,fiber in enumerate(Fibers):
        fiber.spectrum = norm[i,:]            

def get_indices(image, Fibers, fsize, cols=None):
    ygrid, xgrid = np.indices(image.shape)                                           
    a,b = image.shape
    diff = np.zeros((len(Fibers),))
    if cols is None:
        cols = np.arange(int(b/10.),int(9.*b/10.))
    for i,fiber in enumerate(Fibers):
        diff[i] = fiber.trace[cols].max() - fiber.trace[cols].min()
    if Fibers:
        cut_size = int(diff.max() + fsize*2 + 1)
    # Calculate initial fiber model and trace-y reference frame
    for i,fiber in enumerate(Fibers):
        my1 = int(np.max([0,fiber.trace.min()-fsize-1]))
        my2 = int(np.min([a,fiber.trace.max()+fsize+1]))
        
        try:
            ix = ygrid[my1:my2,:] - fiber.trace*np.ones((my2-my1,1))
        except:
            print(my1,my2,cut_size)
            sys.exit(1)
        y,x = np.where(np.abs(ix)<=fsize)
        ind = np.argsort(x)
        x=x[ind]
        y=y[ind]
        fiber.yoff = ix[y,x]
        fiber.yind = y + my1
        fiber.xind = x 

def measure_background(image, Fibers, width=30, niter=3, order=3):
    t = []
    a,b = image.shape
    ygrid,xgrid = np.indices(image.shape)
    ygrid = 1. * ygrid.ravel() / a
    xgrid = 1. * xgrid.ravel() / b
    image = image.ravel()
    s = np.arange(a*b)
    for fiber in Fibers:
        t.append(fiber.D*fiber.yind + fiber.xind)
    t = np.hstack(t)
    t = np.array(t, dtype=int)
    ind = np.setdiff1d(s,t)
    back = biweight_location(image[ind])
    return back
    mask = np.zeros((a*b))
    mask[ind] = 1.
    mask[ind] = 1.-is_outlier(image[ind])
    sel = np.where(mask==1.)[0]
    for i in xrange(niter):
        V = polyvander2d(xgrid[sel],ygrid[sel],[order,order])
        sol = np.linalg.lstsq(V, image[sel])[0]
        vals = np.dot(V,sol) - image[sel]
        sel = sel[~is_outlier(vals)]
    V = polyvander2d(xgrid,ygrid,[order,order])
    back = np.dot(V, sol).reshape(a,b)    
    return back
        
            
def fast_nonparametric_fibermodel(image, Fibers, fsize, bins, sigma, power, 
                                  c1=0.001, c2=0.002, break_fix=4.5, 
                                  col_group=40, fib_group=4, cols=None, 
                                  mask=None, kind='linear', niter=1,
                                  outfile='temp.png', do_fit=True):                                             
    if mask is None:
        mask = np.zeros(image.shape)
    ygrid, xgrid = np.indices(image.shape)                                           
    a,b = image.shape
    
    # Re-normalization so the total model (core+plaw) = 1 for a given column
    I = get_norm_for_gauss_power()
    #pl_norm = get_norm_for_power_law(c1, sigma, a)
    
    # 1st order correction for modelling amplifier instead of spectrograph
    #pl_norm1 = get_norm_for_power_law(c1, sigma, fsize*2)
    #missing_flux = pl_norm - pl_norm1
        
    # Calculate cut_size                                      
    if cols is None:
        cols = np.arange(b)

    get_indices(image, Fibers, fsize)
    for i, fiber in enumerate(Fibers):
        fiber.core = gauss_power(fiber.yoff, sigma, power, 
                                 norm=1./I(sigma, power))
    # Polynomial model in column direction with "missing" flux as total    
    # Calculate spectrum    
    col_grid = cols[::col_group]
    if col_grid[-1] != Fibers[0].D:
        col_grid = np.append(col_grid, Fibers[0].D)
    interp_list = get_interp_list(fsize, bins, kind)
    fib_grid = np.arange(fib_group,len(Fibers)+1,fib_group*2)
    for k in xrange(niter):
        new_fast_norm(image, Fibers, cols, mask)
        xbink=[]
        solk=[]
        xloc=[]
        yloc=[]
        for j,col in enumerate(col_grid[:-1]):
            for i in fib_grid:
                if do_fit:
                    xbin, sol = remeasure_core(image, Fibers, i, fib_group, [col_grid[j],col_grid[j+1]], 
                                               bins, cols, fsize, c1, c2, break_fix, 
                                               interp_list)
                else:
                    xbin = np.linspace(-fsize,fsize, bins)
                    sol = gauss_power(xbin, sigma, power, norm=1./I(sigma,power))
                xbink.append(xbin)
                solk.append(sol)
                xloc.append(col_grid[j]/2.+col_grid[j+1]/2.)
                yloc.append(i)
        sv = []
        nv = []
        X = np.zeros((len(xloc),2))
        X[:,0] = np.array(xloc)/(1.*Fibers[0].D)
        X[:,1] = np.array(yloc)/(1.*len(Fibers))
        for i in xrange(bins):
            sv.append(LinearNDInterpolator(X, np.array(solk)[:,i]))
            nv.append(NearestNDInterpolator(X, np.array(solk)[:,i]))
        for i,fiber in enumerate(Fibers):
            fiber.core = np.zeros(fiber.xind.shape)
            fiber.fibmodel = np.zeros((bins, fiber.D))
            x =np.zeros((len(fiber.xind),2))
            x[:,0] = fiber.xind/(1.*fiber.D)
            x[:,1] = i/(1.*len(Fibers))*np.ones((len(fiber.xind),))
            x1 = np.zeros((fiber.D,2))
            x1[:,0] = np.arange(fiber.D)/(1.*fiber.D) 
            x1[:,1] = i/(1.*len(Fibers))*np.ones((fiber.D,))
            for j,iv in enumerate(interp_list):
                si = sv[j](x)
                if np.isnan(si).sum():    
                    si[np.isnan(si)] = nv[j](x[np.isnan(si)])
                si1 = sv[j](x1)
                if np.isnan(si1).sum():    
                    si1[np.isnan(si1)] = nv[j](x1[np.isnan(si1)])    
                fiber.core += iv(fiber.yoff)*si
                fiber.fibmodel[j,:] = 1.*si1               
    new_fast_norm(image, Fibers, np.arange(b), mask)
    
    for i,fiber in enumerate(Fibers):
        subtract_neighbor_fibers(image, Fibers, i, 1)
        fiber.profile = fiber.subtracted / fiber.spectrum[fiber.xind]
     
    
def fit_fibermodel_nonparametric_bins(image, xgrid, ygrid, Fibers, PL, fib=0, 
                                      xlow=0, xhigh=1032, plot=False, 
                                      group=4, bins=11, niter=3, debug=False,
                                      outfolder=None, sol=None, binx=None,
                                      do_fit=True):
    '''
    This is the workhorse function for fibermodel fitting.  It fits
    a single empirical fiber model for a group of fibers and columns.  
    
    :param image:
        Amplifier image
    :param xgrid:
        Column indices
    :param ygrid:
        Row indices
    :param Fibers:
        List of fiber class object for each fiber
    :param fib:
        Fiber number for fitting
    :param xlow:
        The low column value to use when fitting
    :param xhigh:
        The high column value to use when fitting
    :param plot:
        Plot the fit
    :param fsize:
        Region in which fiber model is defined and at the end set to zero.
        The region is [-fsize, fsize] in pixels.  
    :param group:
        Group size of fibers
    :param bins:
        The initial number of bins used to describe the 
        empirical/non-parametric profile.  Two other bins are added near
        the peak for common asymetries that are in the data but not in the
        initial_profile that defines the x-values of the bins.  This assures
        more accurate modeling of the peak of the profiles.
    :param niter:
        Number of iterations between normalization and fiber model.  Each are 
        fit independently with the other fixed.
    :param debug:
        General timing and debugging
    :param outfolder:
        File folder location of the plot
    :param sol:
        Initial solution for the fiber model to start with
    :param binx:
        Location of the bins in the fiber direction 
        
    '''
    buf = 0.
    if debug:
        t1 = time.time()    
    # Get initial model        

    # Normalize the profiles to 1.
    sol /= np.sum(sol[:-1] * np.diff(binx))
    
    # The lowest and highest fiber to be used in modeling
    lowfib = np.max([0,fib-group/2])
    highfib = np.min([len(Fibers)-1,fib+group/2])
    
    # get low y and high y for fit
    mn1 = np.min(Fibers[lowfib].trace[xlow:xhigh])-buf
    mx1 = np.max(Fibers[lowfib].trace[xlow:xhigh])+buf
    mn2 = np.min(Fibers[highfib].trace[xlow:xhigh])-buf
    mx2 = np.max(Fibers[highfib].trace[xlow:xhigh])+buf
    ylow = int(np.min([mn1,mn2,mx1,mx2]))
    yhigh = int(np.max([mn1,mn2,mx1,mx2]))+1
    ymean = np.mean(Fibers[fib].trace[xlow:xhigh])
    # Define cutouts and fibers
    fibers = Fibers[lowfib:highfib+1]
    x = xgrid[ylow:yhigh,xlow:xhigh].ravel()
    y = ygrid[ylow:yhigh,xlow:xhigh].ravel()
    z = image[ylow:yhigh,xlow:xhigh].ravel()
    P = PL[ylow:yhigh,xlow:xhigh,lowfib:highfib+1].reshape((yhigh-ylow)*(xhigh-xlow),(highfib+1-lowfib))
    # Dummy error
    #zerr = np.sqrt(image[ylow:yhigh,xlow:xhigh].ravel())
    
    # selection within cutout
    ycutl = ((Fibers[lowfib].trace[xlow:xhigh]-buf)*np.ones((yhigh-ylow,1))).ravel()
    ycuth = ((Fibers[highfib].trace[xlow:xhigh]+buf)*np.ones((yhigh-ylow,1))).ravel()
    sel = np.where((y>=ycutl) * (y<=ycuth))[0]    
    # Create empty arrays for the fibermodel weights in each given pixel from 
    Fl = np.zeros((len(y), bins, len(fibers)))
    fun = np.zeros((bins,))
    for i,fiber in enumerate(fibers):
        ytrace = (fiber.trace[xlow:xhigh]*np.ones((yhigh-ylow,1))).ravel()
        ix = y-ytrace
        for j in xrange(bins):
            fun[j] = 1.0
            Fl[:,j,i] = np.interp(ix,binx,fun,left=0.0,right=0.0)
            fun[j] = 0.0
        
    F = Fl.sum(axis=2)
    
    # Solve for fiber model iteratively with the normalization
    for i in xrange(niter):
        init_model = np.zeros((len(x),len(fibers)))
        full_model = np.zeros((len(x),bins))
        plaw_model = np.zeros((len(x)))
        flat = np.zeros((len(x),))
        #flat_err = np.zeros((len(x),))
        norm = np.zeros((len(fibers),xhigh-xlow))
        normfits = np.zeros((len(x),len(fibers)))
        for j,v in enumerate(np.arange(xlow,xhigh)):
            xsel = np.intersect1d(np.where(x==v)[0],sel)
            for k,fiber in enumerate(fibers):
                init_model[xsel,k] = np.dot(Fl[xsel,:,k],sol) + P[xsel,k]

            # Solve for the normalization of the number of fibers  
            try:          
                norm[:,j] = lstsq(init_model[xsel,:],z[xsel])[0]
            except ValueError:
                print(lowfib, highfib, ylow, yhigh, xlow, xhigh, ycutl, ycuth)
                sys.exit(1)
            normfits[xsel,:] = np.ones((len(xsel),1)) * norm[:,j]
            for k in xrange(bins):
                full_model[xsel,k] = (Fl[xsel,k,:] 
                                      * normfits[xsel,:]).sum(axis=1)
            plaw_model[xsel] = (P[xsel,:] * normfits[xsel,:]).sum(axis=1)

        
        try:
            if do_fit:
                sol1 = nnls(full_model[sel,1:-1], z[sel] - plaw_model[sel])[0]
            else:
                sol1 = sol[1:-1]
        except ValueError:
            print(full_model[sel])
            print(plaw_model[sel])
            print(z[sel])
            print("nnls for fiber model solution failed on iteration %i" %i)
            sys.exit(1)
        sol = np.hstack([0.0,sol1,0.0])
        sol /= np.sum(sol[:-1] * np.diff(binx))
    
    if debug:
        t2 = time.time()
        print("Solution took: %0.3f s" %(t2-t1))        
    if plot:
        
        model = np.dot(F,sol) + P.sum(axis=1)
        Fp = np.zeros((len(x),len(fibers)))
        for i in xrange(len(fibers)):
            Fp[:,i] = np.dot(Fl[:,:,i],sol) + P[:,i]
        W = (normfits * Fp).sum(axis=1) / Fp.sum(axis=1)
        flat = np.where(W!=0, z / W, 0.0)
        fig = plt.figure(figsize=(8,6))
        fibplot = plt.axes([0.1,0.55,0.8,0.4])
        implot = plt.axes([0.1,0.1,0.28,0.4])
        normplot = plt.axes([.55,.1,0.35,0.4])
        ytrace = (Fibers[fib].trace[xlow:xhigh]*np.ones((yhigh-ylow,1))).ravel()
        fibplot.plot([ylow-ymean,yhigh-ymean],[0,0],'k-',lw=3,color=[.6,.3,.6],zorder=1)
        fibplot.scatter(y-ytrace, flat, c=[1.0,0.31,0.31], edgecolor='none', 
                        alpha=0.5, s=8)
        for i in xrange(len(Fl[0,0,:])):
            lmodel = np.dot(Fl[:,:,i],sol) + P[:,i]
            fibplot.scatter(y-ytrace, lmodel, c=[0.41,0.41,0.91],edgecolor='none',
                    alpha=0.5,s=8)
        fibplot.scatter(y-ytrace, model, c=[0.11,0.11,0.91],edgecolor='none',
                    alpha=0.5,s=5)
        fibplot.scatter(y-ytrace, flat - model, c=[0.41,0.41,
                                                   0.41],edgecolor='none',
                    alpha=0.5,s=5, zorder=10)
        fibplot.axis([ylow-ymean, yhigh-ymean, -0.01, 0.21])
        cmap = plt.get_cmap('RdBu_r')
        
        cax = implot.scatter(x[sel], y[sel], c=flat[sel]/model[sel], 
                             cmap=cmap, edgecolor='none', s=5, vmin=0.95, 
                             vmax=1.05,marker='s')
        for i in xrange(lowfib,highfib+1):
            yt = Fibers[i].trace[xlow:xhigh]
            implot.plot(np.arange(xlow,xhigh), yt, c=[0.3,0.8,0.32], alpha=0.4)
        implot.axis([xlow, xhigh, ylow, yhigh])
        implot.axis('scaled')
        cbaxes = fig.add_axes([0.4, 0.1, 0.03, 0.4]) 
        cb = plt.colorbar(cax, cax = cbaxes)            
        cmap = matplotlib.cm.get_cmap('viridis')
        mx = highfib+1 - lowfib
        for i in xrange(mx):
            normplot.plot(np.arange(xlow,xhigh), norm[i,:], 
                          c=cmap((1.*i)/(mx-1)))
        normplot.set_xlim([xlow,xhigh])
        fig.savefig(op.join(outfolder, "fit_bins_%s_%03d_%04d_%04d.png" 
                    %(Fibers[fib].basename, fib, xlow, xhigh)), dpi=150)
        plt.close(fig)    
    return sol


def get_norm_nonparametric_fast(image, Fibers, cols=None, mask=None,
                                plaw_coeff1=0.0004):
    plaw_coeff = np.array([plaw_coeff1, 0.5, 0.15, 1.0])

    bins=len(Fibers[0].binx)
    a,b = image.shape
    if mask is None:
        mask = np.zeros(image.shape)
    ygrid,xgrid = np.indices(image.shape)
    fun = np.zeros((bins,))
    y=ygrid[:,0]
    Fl = np.zeros((len(y), bins))
    P = build_powerlaw(ygrid, Fibers, plaw_coeff, cols)
    init_model = np.zeros((len(y),len(Fibers)))
    if cols is None:
        cols = np.arange(b)
    norm = np.zeros((len(Fibers),b))
    for col in cols:
        for i,fiber in enumerate(Fibers):
            ix = y-fiber.trace[col]
            for j in xrange(bins):
                fun[j] = 1.0
                Fl[:,j] = np.interp(ix,fiber.binx,fun,left=0.0,right=0.0)
                fun[j] = 0.0
            init_model[:,i] = np.dot(Fl, fiber.fibmodel[col,:])+P[:,col,i]
        norm[:,col] = lstsq(init_model[mask[:,col]==0,:], 
                            image[mask[:,col]==0,col])[0]
    return norm
 

def get_model_image(image, fibers, prop, debug=False):
    '''
    Produce an amplifier image of a given fiber property, prop.  For example if 
    prop==spectrum, then this function produces a model of the image.
    If instead prop==sky_spectrum, than this function produces an image
    of the sky.  
    
    :param image:
        Amplifier image
    :param fibers:
        List of fiber class object for each fiber
    :param prop:
        A property of the fiber class to be used for the model image.
    :param debug:
        Timing and debugging
    
    '''
    model = np.zeros(image.shape)
    for i,fiber in enumerate(fibers):
        yind = fiber.yind
        xind = fiber.xind
        model[yind,xind] += (getattr(fiber,prop)[xind] * 1.00
                               * fiber.core)

    return model   
    
def check_fiber_trace(image, Fibers, outfile, xwidth=75., ywidth=75.):
    '''
    Plot of the amplifier image overlayed with the trace for 
    the top/middle/bottom in the fiber and wavelength direction (3 x 3)
    
    :param image:
        Amplifier image
    :param Fibers:
        List of fiber class object for each fiber
    :param outfile:
        File name to be saved for this plot
    :param xwidth:
        Width in wavelength/column direction of each box
    :param ywidth:
        Width in fiber/row direction of each box
        
    '''
    ylen,xlen = image.shape
    ypos, xpos = np.indices((ylen,xlen))
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    # initial plot position    
    fig = plt.figure(figsize=(12, 12))    
    pos = 0
    # For plotting purposes of the WaveCheck_{SPECID}_{SIDE}.pdf plots
    plots = np.arange(1,10)
    cmap = plt.get_cmap('Blues')
    for i in [1, 0.5, 0]:
        for j in [0.1, 0.5, 0.9]:
            sub = fig.add_subplot(3, 3, plots[pos])
            pos += 1
            minx = int(j * -1 * xwidth + j * (xlen - 1.))
            maxx = int((j-1) * -1 * xwidth + j * (xlen - 1.))
            miny = int(i * -1 * ywidth + i * (ylen - 1.))
            maxy = int((i-1) * -1 * ywidth + i * (ylen - 1.))
            zarr = image[miny:maxy,minx:maxx].ravel()
            good = np.where(is_outlier(zarr)<1)[0]
            vmin = np.min(zarr[good])
            vmax = np.max(zarr[good])
            sub.imshow(image[miny:maxy,minx:maxx],origin='lower',
                       extent=[minx,maxx,miny,maxy],interpolation='nearest',
                       cmap=cmap, vmin=vmin, vmax=vmax)
            for k in xrange(len(Fibers)):
                sub.plot(np.arange(xlen), Fibers[k].trace+0.5, color=[1.0,0.4,0.35],
                         linewidth=2)
            plt.axis([minx, maxx, miny, maxy])
    fig.savefig(outfile)
    plt.close(fig)     

def check_fiber_profile(image, Fibers, outfile, fsize,
                        xwidth=50., yplot_high=0.25, 
                        yplot_low=-0.05, fib_group=4):
    '''
    Plot of the fiber profiles for the top/middle/bottom in the fiber and 
    wavelength direction (3 x 3)
    
    :param image:
        Amplifier image
    :param Fibers:
        List of fiber class object for each fiber
    :param outfile:
        File name to be saved for this plot
    :param fiber_sel:
        Fibers to plotted
    :param xwidth:
        Width in wavelength/column direction of each box
    :param ywidth:
        Width in fiber/row direction of each box
    :param yplot_high:
        Ylim in the plot (related to the normalized fiber profile height)
    :param yplot_low:
        Ylim in the plot
    :param plotbuf:
        Buffer in fiber direction for plotting
        
    '''
    mta = biweight_location(np.diff([fiber.trace[int(fiber.D/2)] 
                                     for fiber in Fibers]))
    fig = plt.figure(figsize=(12, 12)) 
    a,b = image.shape
    pos = 0
    plots = np.arange(1,10)
    for i in [0.9, 0.5, 0.1]:
        for j in [0.2, 0.5, 0.8]:
            sub = fig.add_subplot(3, 3, plots[pos])
            pos += 1
            minx = (j*b-xwidth/2)
            maxx = (j*b+xwidth/2)
            x = []
            y = []
            mt = []
            mxk = []
            for fib in np.arange(int(i*len(Fibers))-fib_group,
                                 int(i*len(Fibers))+fib_group+1):
                mt.append(Fibers[fib].trace[int(maxx/2.+minx/2.)])
                sel = np.where((Fibers[fib].xind>minx)
                                *(Fibers[fib].xind<maxx))[0]
                x.append(Fibers[fib].yoff[sel])
                y.append(Fibers[fib].core[sel])
                
                sub.scatter(Fibers[fib].yoff[sel],
                            Fibers[fib].profile[sel], alpha=0.03, 
                            edgecolor='none',color=[1.0,0.4,0.4],s=70)
                sub.scatter(Fibers[fib].yoff[sel],
                            Fibers[fib].profile[sel]
                            -Fibers[fib].core[sel],
                            alpha=0.03, edgecolor='none', 
                            color=[0.4, 0.4, 0.4])
                sub.scatter(Fibers[fib].yoff[sel],
                            Fibers[fib].core[sel],
                            color=[0.2, 0.2, 0.8], edgecolor='none',
                            alpha=0.2, s=20)
                mxk.append(Fibers[fib].yoff[sel].max())
            x = np.hstack(x)
            y = np.hstack(y)
            if len(mt)>3:
                mtd = biweight_location(np.diff(mt))
            else:
                mtd = mta*1.
            minsel = np.argmin(np.abs(x-mtd/2.))
            maxsel = np.argmin(np.abs(x))
            peak = y[maxsel]
            trough = y[minsel]
            contrast = (peak - trough*2.) / (peak + trough*2.)
            sub.text(-fsize+1, yplot_high-0.03, 
                     "Fiber: %03d, X: %04d-%04d" %(int(i*len(Fibers)),
                                                   minx+1,maxx+1), 
                     fontsize=14)
            sub.text(-fsize+1, yplot_high-0.06, 
                     "Contrast: %0.2f" %(contrast), 
                     fontsize=14)
            sub.set_xlim([-fsize, fsize])
            sub.set_ylim([yplot_low, yplot_high])                
            sub.plot([-fsize,fsize],[0,0],'r-',lw=3)
    fig.savefig(outfile)
    plt.close(fig)
    
def check_wavelength_fit(Fibers, sun, outfile, fiber_sel=[107,58,5], 
                        xwidth=125., fwidth=10, smooth_length=21):
    '''
    Plot of the wavelength solution for the top/middle/bottom in the fiber and 
    wavelength direction (3 x 3)
    
    :param Fibers:
        List of fiber class object for each fiber
    :param sun:
        Solar template spetrum (two column array)
    :param outfile:
        File name to be saved for this plot
    :param fiber_sel:
        Fibers to plotted
    :param xwidth:
        Width in wavelength/column direction of each box
    :param fwidth:
        Number of fibers included on each side of the fiber selected
    :param smooth_length:
        Biweight filter smoothing length for normalizing the fiber spectra.
        This could be inconsistent with the template and the fitting.
        Be careful.
        
    '''
    xlen = len(Fibers[0].wavelength)
    flen = len(Fibers)
    fig = plt.figure(figsize=(12, 12))    
    pos = 0
    plots = np.arange(1,10)                       
    for i in fiber_sel:
        for j in [0.15, 0.5, .85]:
            sub = fig.add_subplot(3, 3, plots[pos])
            pos += 1
            minx = np.max([0, int(j * -1 * xwidth + j * (xlen - 1.))])
            maxx = np.min([xlen, int((j-1) * -1 * xwidth + j * (xlen - 1.))])
            minf = np.max([0, i-fwidth])
            maxf = np.min([flen, i+fwidth+1])
            minw = Fibers[i].wavelength[minx] - 10.
            maxw = Fibers[i].wavelength[maxx] + 10.
            xl = np.searchsorted(sun[:,0], minw)
            xh = np.searchsorted(sun[:,0], maxw)
            sub.step(sun[xl:xh,0], sun[xl:xh,1],'r-', lw=3, where='mid')
            for fib in xrange(int(minf),int(maxf)):
                y = Fibers[fib].spectrum
                y = biweight_filter(y, smooth_length) / y
                sub.step(Fibers[fib].wavelength[minx:maxx], y[minx:maxx], 'b-',
                         alpha=0.1, where='mid', lw=2)
            sub.set_xlim([minw+10, maxw-10])
            ypl = np.percentile(sun[xl:xh,1],1)
            yph = np.percentile(sun[xl:xh,1],99)
            ran = yph - ypl
            ylow = ypl - 0.2*ran
            yhigh = ypl + 1.2*ran
            sub.text(minw+15,yhigh-0.1*ran, 'Fibers: %03d-%03d' %(minf+1,maxf),
                     fontsize=12)
            sub.set_ylim([ylow, yhigh])
    fig.savefig(outfile)
    plt.close(fig)  


def calculate_significance(Fibers, image, error, oversample_value=3, sigma=1.3, 
                           wave_step=4, sigthresh=4, cosmic_avoidance=4, 
                           interactive=False):
    '''
    Calculates the significance of map from extracted spectra.  
    '''
    # Oversample and Convolve
    for fiber in Fibers:
        s = np.diff(fiber.wavelength) / oversample_value
        t = np.hstack([s, s[-1]])
        wave_list = []
        for i in np.arange(oversample_value):
            wave_list.append([fiber.wavelength + i*t])
        fiber.oversampled_wave = np.sort(np.hstack(wave_list)).ravel()
        oversampled_spec = np.interp(fiber.oversampled_wave, fiber.wavelength, 
                                     fiber.sky_subtracted).ravel()
        fiber.oversampled_ftf = np.interp(fiber.oversampled_wave, fiber.wavelength, 
                                     fiber.fiber_to_fiber).ravel()
        kernel = Gaussian1DKernel(stddev=sigma*oversample_value)
        fiber.convolved_spec = convolve(oversampled_spec, kernel)
    spectrum = np.array([fiber.convolved_spec for fiber in Fibers])
    wavelength = np.array([fiber.oversampled_wave for fiber in Fibers])
    ftf = np.array([fiber.oversampled_ftf for fiber in Fibers])
    trace = np.array([fiber.trace for fiber in Fibers])
    
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
    Sig = spectrum/(V*np.sqrt(ftf))
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
    
    if interactive:
        import pyds9
        ds9 = pyds9.DS9()
        ds9.set_np2arr(image)
        ds9.set('scale mode zscale')
        fibind, xind = np.where(Sig>sigthresh)
        fib_dict = {}
        for i,fib in enumerate(Fibers):
            fib_dict[i] = 0
        for fib, xind in zip(fibind,xind):
            fib_dict[fib] +=1
            if fib_dict[fib] < 50:
                x = 1.*xind/oversample_value + 1
                y = Fibers[fib].trace[int(x)-1]
                s = Sig[fib,xind]
                ds9.set('region','image; circle %f %f %f # color=blue' %(x,y,s))
    return Sig, wavelength
   
def get_wavelength_offsets(fibers, binsize, wvstep=0.2, colsize=300):
    bins = np.hstack([np.arange(0,fibers[0].D,binsize),fibers[0].D])
    mn = np.min([fiber.wavelength.min() for fiber in fibers])
    mx = np.min([fiber.wavelength.max() for fiber in fibers])
    yf = []
    wv = np.arange(mn,mx,wvstep)
    pd = []
    nbins = len(bins)-1
    for fiber in fibers:
        y = fiber.spectrum
        w = fiber.wavelength
        x = np.arange(fiber.D)
        s = -999*np.ones((nbins,fiber.D))
        for i in np.arange(nbins):
            x0 = int(np.max([0,(bins[i]+bins[i+1])/2 - colsize/2]))
            x1 = int(np.min([fibers[0].D,(bins[i]+bins[i+1])/2 + colsize/2]))
            
            s[i,x0:x1], flag = polynomial_normalization(x, y, x0, x1)
        s = np.ma.array(s,mask=s==-999.)
        ym = s.mean(axis=0)
        ym = biweight_filter(ym, 5, ignore_central=1)
        pr, ph = find_maxima(w, -1.*(y/ym-1)+1, interp_window=3)
        pd.append(pr[~np.isnan(pr)])
        yf.append(np.interp(wv,w,-1.*(y/ym-1)+1,left=0.0,right=0.0))
    master = biweight_location(yf,axis=(0,))
    pr,ph = find_maxima(wv, master, y_window=10, interp_window=2.5, 
                        repeat_length=3, first_order_iter=5, 
                        max_to_min_dist=15)
    sel = (ph>1.03) * (~np.isnan(pr))
    pr = pr[sel]
    ph = ph[sel]
    wk = []
    for p,fiber in zip(pd,fibers):
        
        c = p[:,None] - pr
        a = np.where(np.abs(c)<4.)
        x = pr[a[1]]
        y = p[a[0]] - pr[a[1]]
        sel = np.argsort(x)
        x = x[sel]
        y = y[sel]
        m = medfilt(y,15)
        good = ~is_outlier(y-m) 
        p0 = np.polyfit(x[good],y[good],3)
        wk.append(np.polyval(p0,fiber.wavelength))
    return wk
            
def polynomial_normalization(x, y, x0, x1, lowsig=0.5, highsig=3.,
                                       niter=5, order=3):
    xf = x[x0:x1]
    yf = y[x0:x1]
    p0 = np.polyfit(xf,yf,1)
    ym = np.polyval(p0,xf)
    diff = ym - yf
    mad = np.median(np.abs(diff))
    mzs = 0.6745 * diff / mad
    sel = (mzs<lowsig)*(mzs>-highsig)
    for i in np.arange(niter):
        p0 = np.polyfit(xf[sel],yf[sel],1)
        ym = np.polyval(p0,xf)
        diff = ym - yf
        mad = np.median(np.abs(diff))
        mzs = 0.6745 * diff / mad
        sel = (mzs<lowsig)*(mzs>-highsig)
    p0 = np.polyfit(xf[sel],yf[sel],order)
    return np.polyval(p0,xf), sel
    
def bspline_x0(x, y=None, nknots=21, nknots2=21, p=3):
    if y is not None:
        s = np.cumsum(np.abs(np.diff(np.diff(y))))
        s = np.hstack([0,s/s.max(),1])
        v = np.interp(np.linspace(0,1,nknots), s, 
                      (x-x.min())/(x.max()-x.min()))
        v = np.sort(np.hstack([np.linspace(0,1,nknots2),
                               v]))
    else:
        v = np.linspace(0,1,nknots)
    k = splinelab.augknt(v, p)  # add endpoint repeats as appropriate for spline order p
    B = Bspline(k, p)
    xi =  (x-x.min())/(x.max()-x.min()+0.1)
    c = np.array([B(xp) for xp in xi])
    return B, c

def fit_continuum_sky(wv, sky, fil_len=95, func=np.array):
    skym_s = 1.*sky
    sky_sm = savgol_filter(skym_s, fil_len, 1)
    for i in np.arange(5):
        mad = np.median(np.abs(sky - sky_sm))
        outlier = func(sky - sky_sm) > 1.5 * mad
        skym_s = 1.*sky
        skym_s[outlier] = np.interp(wv[outlier], wv[~outlier],
                                    sky_sm[~outlier])
        sky_sm = savgol_filter(skym_s, fil_len, 1)
    return sky_sm
