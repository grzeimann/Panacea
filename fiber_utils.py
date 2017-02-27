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
from scipy.linalg import lstsq
import matplotlib.pyplot as plt
import matplotlib
import os.path as op
import numpy as np
import time
import sys

plaw_coeff = np.array([0.0004,0.5,0.15,1.0])


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
                    fixscale=False): 
    '''
    Fitting function for linear wavelength solution of a small column section
    of one fiber.  
    '''
    if fixscale:
        def f(params, sel1):
            wv = init_scale * np.arange(D) + params[0]
            model = np.interp(wv[xi:xe+1], sun_wave, ysun, left=0.0, right=0.0)
            return model[sel1] - data[sel1]
        params0 = np.array([init_wave0])
        sel = is_outlier(f(params0,np.ones(data.shape,dtype=bool)))<1
        sol = scipy.optimize.leastsq(f, params0, args=(sel))[0]
        chi2 = f(sol, sel+True )**2
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


def calculate_wavelength_chi2(x, y, solar_spec, smooth_length=21,
                         init_lims=None, order=3, init_sol=None, debug=False,
                         interactive=False, nbins=21, wavebuff=100, 
                         plotbuff=70, fixscale=False):
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
    y = biweight_filter(y, smooth_length) / y
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
            if init_sol is None:
                xwv0 = np.arange(wv0-20,wv0+20,5)
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
                                  fixscale=fixscale)
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
                ax.set_ylim([-.1*rn + mn, 1.1*rn+mn])
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

def get_trace_from_image(image, y_window=3, x_window=5, repeat_length=2,
                         order=3, mx_cut=0.1, max_to_min_dist=5., debug=False,
                         first_order_iter=5, interp_window=3.):
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
    xc = np.arange(b)[x_window:(b-1-x_window)]
    if debug:
        t1 = time.time()
    for i in xc:
        # Average 2*x_windwow columns at a time
        y = biweight_location(image[:,(i-x_window):(i+x_window+1)], axis=(1,))
        
        # Set empty arrays for maxima and minima found in 2*y_window searches
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


def plaw(xp, plaw_coeff):
    return plaw_coeff[0] / (plaw_coeff[1] + plaw_coeff[2]
              * np.power(abs(xp / 2.5), plaw_coeff[3]))    
    
def fit_fibermodel_nonparametric(image, Fibers, plot=False, fsize=8., 
                                 fiber_group=4, bins=15, col_group=48,
                                 debug=False, use_default=False,
                                 outfolder=None, sigma=2.5, power=2.5):
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
    a,b = image.shape 
    ygrid,xgrid = np.indices(image.shape)                       
    nfibs = len(Fibers) 
    ncols = b / col_group
    so = np.zeros((nfibs, ncols, bins+2))
    xcols = np.arange(col_group/2, int((ncols-1/2.)*col_group)+1, col_group)
    bins1, binx, sol = init_fibermodel(fsize=fsize, bins=bins, sigma=sigma,
                                       power=power)
    PL = np.zeros((a,b,nfibs))
    for i in xrange(b):
        for j,fiber in enumerate(Fibers):
            ix = ygrid[:,i] - fiber.trace[i]
            PL[:,i,j] = plaw(ix, plaw_coeff)
            
            
    if debug:
        t1 = time.time()
    for j in xrange(nfibs):
        for i in xrange(ncols):
            if not use_default:
                sol = fit_fibermodel_nonparametric_bins(image, xgrid, ygrid, 
                                                    Fibers, PL, fib=j, 
                                                    debug=debug, 
                                                    group=fiber_group, 
                                                    bins=bins1,
                                                    xlow=i*col_group, 
                                                    xhigh=(i+1)*col_group, 
                                                    plot=plot,
                                                    outfolder=outfolder,
                                                    sol=sol, binx=binx)                
            so[j,i,:] = sol         
    if debug:
        t2 = time.time()
        print("Fibermodel Solution took: %0.3f s" %(t2-t1))
    if False:#plot:
        for i in xrange(len(sol)):
            mn = np.percentile(so[:,:,i],10)
            mx = np.percentile(so[:,:,i],90)
            plt.figure(figsize=(8,6))
            plt.imshow(so[:,:,i], vmin=mn, vmax=mx, interpolation='nearest', 
                       origin='lower',extent=[0,1032,0,2064])
            plt.colorbar()
    return so, xcols, binx

def get_norm_nonparametric(image, Fibers, fiber_group=4, 
                           debug=False, mask=None):
    '''
    This builds the normalization (aka spectrum) for each fiber and mimicks
    the fibermodel function above but instead of fitting the profile, uses
    an already fitted profile and solves for the spectrum.
    
    :param image:
        Amplifier image
    :param Fibers:
        List of fiber class object for each fiber
    :param fiber_group:
        Because the grid is built with fiber_group, this is an input just as
        it was for fibermodel but allows the normalization to be solved for
        over a different grid if someone wants to use more/less fibers
    :param debug:
        All-purpose boolean for debugging
    :param mask:
        Image with -1 as a mask for ignoring pixels when solving for the 
        spectrum.  It is not clear yet if a whole fiber is flagged over a large
        number of columns, whether a reasonable value or will lead to nonsense.
    '''
    a,b = image.shape 
    if mask is None:
        mask = np.zeros(image.shape)
    ygrid,xgrid = np.indices(image.shape)                       
    nfibs = len(Fibers) 
    ncols = b 
    norm = np.zeros((nfibs, ncols))
    if debug:
        t1 = time.time()
    for j in xrange(nfibs):
        norm[j,:] = get_norm_nonparametric_bins(image, mask, xgrid, ygrid, 
                                                    Fibers, fib=j, 
                                                    group=fiber_group, xlow=0, 
                                                    xhigh=ncols, debug=debug)                
    if debug:
        t2 = time.time()
        print("Solution took: %0.3f s" %(t2-t1))

    return norm   
    
    
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


def get_norm_nonparametric_fast(image, Fibers, cols, mask=None):
    bins=len(Fibers[0].binx)
    a,b = image.shape
    if mask is None:
        mask = np.zeros(image.shape)
    ygrid,xgrid = np.indices(image.shape)
    fun = np.zeros((bins,))
    y=ygrid[:,0]
    Fl = np.zeros((len(y), bins))
    Pl = np.zeros((len(y),))
    init_model = np.zeros((len(y),len(Fibers)))
    norm = np.zeros((len(Fibers),b))
    for col in cols:
        for i,fiber in enumerate(Fibers):
            ix = y-fiber.trace[col]
            for j in xrange(bins):
                fun[j] = 1.0
                Fl[:,j] = np.interp(ix,fiber.binx,fun,left=0.0,right=0.0)
                fun[j] = 0.0
            Pl = plaw(ix, plaw_coeff)
            init_model[:,i] = np.dot(Fl, fiber.fibmodel[col,:])+Pl
        norm[:,col] = lstsq(init_model[mask[:,col]==0,:], 
                            image[mask[:,col]==0,col])[0]
    return norm
        
def get_norm_nonparametric_bins(image, mask, xgrid, ygrid, Fibers, fib=0, 
                                      xlow=0, xhigh=1032, fsize=8., 
                                      group=4, debug=False):
    '''
    This is the workhorse function for producing the spectrum of a
    fiber.  It is similar to fit_fibermodel_nonparametric_bins.
    
    :param image:
        Amplifier image
    :param mask:
        Pixels to ignore in fitting for the normalization.  Mask==-1 is ignored
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
    :param debug:
        General timing and debugging

    '''
    if debug:
        t1 = time.time()    
    
    # The lowest and highest fiber to be used in modeling
    lowfib = np.max([0,fib-group/2])
    highfib = np.min([len(Fibers)-1,fib+group/2])
    
    # get low y and high y for fit
    mn1 = np.min(Fibers[lowfib].trace[xlow:xhigh])
    mx1 = np.max(Fibers[lowfib].trace[xlow:xhigh])
    mn2 = np.min(Fibers[highfib].trace[xlow:xhigh])
    mx2 = np.max(Fibers[highfib].trace[xlow:xhigh])
    ylow = int(np.min([mn1,mn2,mx1,mx2]))
    yhigh = int(np.max([mn1,mn2,mx1,mx2]))+1
    # Define cutouts and fibers
    fibers = Fibers[lowfib:highfib+1]
    fid = np.where(np.arange(lowfib,highfib+1)==fib)[0]
    x = xgrid[ylow:yhigh,xlow:xhigh].ravel()
    y = ygrid[ylow:yhigh,xlow:xhigh].ravel()
    z = image[ylow:yhigh,xlow:xhigh].ravel()
    msk = mask[ylow:yhigh,xlow:xhigh].ravel()
    
    # selection within cutout
    ycutl = (Fibers[lowfib].trace[xlow:xhigh]*np.ones((yhigh-ylow,1))).ravel()
    ycuth = (Fibers[highfib].trace[xlow:xhigh]*np.ones((yhigh-ylow,1))).ravel()
    sel = np.where((y>=ycutl) * (y<=ycuth) * (msk==0))[0]
    # Create empty arrays for the fibermodel weights in each given pixel from
    dummy, bins = fibers[0].fibmodel.shape
    binx = fibers[0].binx
    Fl = np.zeros((len(y), bins, len(fibers)))
    Pl = np.zeros((len(y),len(fibers)))
    fun = np.zeros((bins,))
    i = 0
    for fiber in fibers:
        ytrace = (fiber.trace[xlow:xhigh]*np.ones((yhigh-ylow,1))).ravel()
        ix = y-ytrace
        for j in xrange(bins):
            fun[j] = 1.0
            Fl[:,j,i] = np.interp(ix,binx,fun,left=0.0,right=0.0)
            fun[j] = 0.0
        Pl[:,i] = plaw(ix, plaw_coeff)
        i+=1
        
    
    # Solve for the normalization
    init_model = np.zeros((len(x),len(fibers)))
    norm = np.zeros((len(fibers),xhigh-xlow))
    for j,v in enumerate(np.arange(xlow,xhigh)):
        xsel = np.intersect1d(np.where(x==v)[0],sel)
        if len(xsel):
            k = 0
            for fiber in fibers:
                init_model[xsel,k] = (np.dot(Fl[xsel,:,k],fiber.fibmodel[j,:]) 
                                      + Pl[xsel,k])
                k+=1
            # Solve for the normalization of the number of fibers
            norm[:,j] = lstsq(init_model[xsel,:],z[xsel])[0]
        else:
            norm[:,j] = 0.0
    if debug:
        t2 = time.time()
        print("Fiberextract solution for Fiber %i took: %0.3f s" %(fib, t2-t1))  
    return norm[fid,:]
 

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
    if debug:
        t1 = time.time()    
    
    # Create empty arrays for the fibermodel weights in each given pixel from
    dummy, bins = fibers[0].fibmodel.shape
    binx = fibers[0].binx
    low = binx.min()-1
    high = binx.max()+1
    fun = np.zeros((bins,))
    def plaw(xp, plaw_coeff):
        return plaw_coeff[0] / (plaw_coeff[1] + plaw_coeff[2]
                  * np.power(abs(xp), plaw_coeff[3]))
    model = np.zeros(image.shape)
    a,b = image.shape
    y = np.arange(a,dtype=float)
    for i in xrange(b):
        Fl = np.zeros((a, bins))
        for fib, fiber in enumerate(fibers):
            ix = y - fiber.trace[i]
            li = np.searchsorted(ix,low)
            hi = np.searchsorted(ix,high)
            for j in xrange(bins):
                fun[j] = 1.0
                Fl[li:hi,j] = np.interp(ix[li:hi],binx,fun,left=0.0,right=0.0)
                fun[j] = 0.0
            model[li:hi,i] += (getattr(fiber,prop)[i] * 1.00
                               * (np.dot(Fl[li:hi,:],fiber.fibmodel[i,:]) 
                                  + plaw(ix[li:hi], plaw_coeff)))
    
    if debug:
        t2 = time.time()
        print("Solution for model image took: %0.3f s" %(t2-t1))  
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
        for j in [0, 0.5, 1]:
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

def check_fiber_profile(image, Fibers, outfile, fiber_sel=[5,58,107], 
                        xwidth=75., ywidth=75., yplot_high=0.25, 
                        yplot_low=-0.01, plotbuf=4, neighbor_fibers=3):
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
    ylen,xlen = image.shape
    ypos, xpos = np.indices((ylen,xlen))
    nfibs = len(Fibers)
    # initial plot position    
    fig = plt.figure(figsize=(12, 12))    
    pos = 0
    # For plotting purposes of the WaveCheck_{SPECID}_{SIDE}.pdf plots
    plots = np.arange(1,10)
    fiber_sel = np.sort(np.array(fiber_sel))[::-1]
    dummy, bins = Fibers[0].fibmodel.shape
    binx = Fibers[0].binx
    low = binx.min()-8
    high = binx.max()+8
    fun = np.zeros((bins,))
    PL = np.zeros((ylen,xlen,nfibs))
    for i in xrange(xlen):
        for j,fiber in enumerate(Fibers):
            ix = ypos[:,i] - fiber.trace[i]
            PL[:,i,j] = plaw(ix, plaw_coeff)
    for i in fiber_sel:
        for j in [0, 0.5, 1]:
            sub = fig.add_subplot(3, 3, plots[pos])
            pos += 1
            minx = int(j * -1 * xwidth + j * (xlen - 1.))
            maxx = int((j-1) * -1 * xwidth + j * (xlen - 1.))
            lowfib = np.max([0,i-neighbor_fibers])
            highfib = np.min([len(Fibers)-1,i+neighbor_fibers])
            fibers = Fibers
            # get low y and high y for fit
            mn1 = np.min(Fibers[lowfib].trace[minx:maxx])-plotbuf
            mx1 = np.max(Fibers[lowfib].trace[minx:maxx])+plotbuf
            mn2 = np.min(Fibers[highfib].trace[minx:maxx])-plotbuf
            mx2 = np.max(Fibers[highfib].trace[minx:maxx])+plotbuf
            miny = np.max([int(np.min([mn1,mn2,mx1,mx2])),0])
            maxy = np.min([ylen, int(np.max([mn1,mn2,mx1,mx2]))+1])
                        
            indx = (maxx-minx) / 2
            yhigh = Fibers[highfib].trace[indx] - Fibers[i].trace[indx]+plotbuf
            ylow = Fibers[lowfib].trace[indx] - Fibers[i].trace[indx]-plotbuf

            model = np.zeros(image.shape)
            flat = np.zeros(image.shape)
            normfits = np.zeros((image.shape + (len(Fibers),)))
            lmodel = np.zeros((image.shape + (len(Fibers),)))
            for k in xpos[0,minx:maxx]:
                Fl = np.zeros((ylen, bins))
                for fib, fiber in enumerate(fibers):
                    ix = ypos[:,k] - fiber.trace[k]
                    li = np.searchsorted(ix,low)
                    hi = np.searchsorted(ix,high)
                    ytrace = ypos[li:hi,k] - Fibers[i].trace[k] 
                    for l in xrange(bins):
                        fun[l] = 1.0
                        Fl[li:hi,l] = np.interp(ix[li:hi], binx, fun, left=0.0, 
                                                right=0.0)
                        fun[l] = 0.0
                    lmodel[li:hi,k,fib] = (np.dot(Fl[li:hi,:],
                                                  fiber.fibmodel[k,:]) 
                                           + PL[li:hi,k,fib])
                    model[li:hi,k] += lmodel[li:hi,k,fib]
                    # TODO proper normalization
                    normfits[li:hi,k,fib] = fiber.spectrum[k] * 1.03
                W = ((normfits[miny:maxy,k,:]*lmodel[miny:maxy,k,:]).sum(axis=1) /
                     lmodel[miny:maxy,k,:].sum(axis=1))
                flat[miny:maxy,k] = np.where(W!=0,image[miny:maxy,k]/W, 0.0)
            ytrace = ypos[miny:maxy,minx:maxx] - Fibers[i].trace[minx:maxx]            
            sub.scatter(ytrace.ravel(), flat[miny:maxy,minx:maxx].ravel(), 
                        c=[1.0,0.31,0.31], edgecolor='none', alpha=0.5, s=8)
            sub.scatter(ytrace.ravel(), model[miny:maxy,minx:maxx].ravel(), 
                         c=[0.11,0.11,0.91],edgecolor='none', alpha=0.5,s=5)
            sub.scatter(ytrace.ravel(), flat[miny:maxy,minx:maxx].ravel() 
                        - model[miny:maxy,minx:maxx].ravel(), 
                        c=[0.41,0.41,0.41], edgecolor='none', alpha=0.5,s=5)
            sub.plot([ylow, yhigh], [0.0, 0.0],'k-',lw=2)
            sub.text(ylow+1, yplot_high-0.03, 
                     "Fiber: %03d, X: %04d-%04d" %(i,minx+1,maxx+1), 
                     fontsize=14)
            sub.set_xlim([ylow, yhigh])
            sub.set_ylim([yplot_low, yplot_high])
    fig.savefig(outfile)
    plt.close(fig) 
    
def check_wavelength_fit(Fibers, sun, outfile, fiber_sel=[5,58,107], 
                        xwidth=50., fwidth=10, smooth_length=21):
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
        for j in [0, 0.5, 1]:
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