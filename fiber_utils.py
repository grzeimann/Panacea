# -*- coding: utf-8 -*-
"""
Fiber Utilities
-----------
To be used in conjuction with IFU reduction code, Panacea


"""

import matplotlib
matplotlib.use('agg')
from utils import biweight_location
from scipy.optimize import nnls
import matplotlib.pyplot as plt
import numpy as np
import time
import sys

def get_trace_from_image(image, y_window=3, x_window=5, repeat_length=2,
                         order=3, mx_cut=0.1, max_to_min_dist=5., debug=False,
                         first_order_iter=5):
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
                xp = np.linspace(lwa-y_window, lwa+y_window, num=50)
                yp = np.interp(xp, x[mn1:mn2], y[mn1:mn2],left=0,right=0)
                lwa = (yp*xp).sum()/yp.sum()
            peaks_refined[j] = lwa
            peaks_height[j] = y[val]
        mh = biweight_location(peaks_height)
        peaks_final = peaks_refined[peaks_height > (mx_cut * mh)]
        allfibers.append(peaks_final) 
    if debug:
        t2 = time.time()
        print("Time Taken for Trace: %0.2f" %(t2-t1))
    return allfibers, xc
    
def fit_fibermodel_nonparametric(image, Fibers, plot=False, fsize=8., 
                                 fiber_group=4, bins=15, col_group=48,
                                 debug=False):
    a,b = image.shape 
    ygrid,xgrid = np.indices(image.shape)                       
    nfibs = len(Fibers) 
    ncols = b / col_group
    so = np.zeros((nfibs, ncols, bins+2))
    xcols = np.arange(col_group/2, int((ncols-1/2.)*col_group)+1, ncols)
    bins, binx, sol = init_fibermodel(fsize=fsize, bins=bins)
    if debug:
        t1 = time.time()
    for j in xrange(nfibs):
        for i in xrange(ncols):
            sol = fit_fibermodel_nonparametric_bins(image, xgrid, ygrid, 
                                                    Fibers, fib=j, 
                                                    group=fiber_group, 
                                                    bins=bins, 
                                                    xlow=i*col_group, 
                                                    xhigh=(i+1)*col_group, 
                                                    plot=plot)
            so[j,i,:] = sol         
    if debug:
        t2 = time.time()
        print("Solution took: %0.3f s" %(t2-t1))
    if plot:
        for i in xrange(len(sol)):
            mn = np.percentile(so[:,:,i],10)
            mx = np.percentile(so[:,:,i],90)
            plt.figure(figsize=(8,6))
            plt.imshow(so[:,:,i], vmin=mn, vmax=mx, interpolation='nearest', 
                       origin='lower',extent=[0,1032,0,2064])
            plt.colorbar()
    return so, xcols, binx
        
    
def init_fibermodel(fsize, bins, sigma=2.7, power=2.5):
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
    binx = np.hstack([binx[:(bins/2)],-.4,binx[bins/2],0.4,binx[((bins/2)+1):]])
    bins+=2
    sol = np.interp(binx,xt,yt,left=0.0,right=0.0)
    sol[0] = 0.0
    sol[-1] = 0.0
    return bins, binx, sol
        
def fit_fibermodel_nonparametric_bins(image, xgrid, ygrid, Fibers, fib=0, 
                                      xlow=0, xhigh=1032, plot=False, fsize=8., 
                                      group=4, bins=11, niter=3, debug=False):
    '''
    : param Fibers:
        list of Fiber class objects (length = number of fibers)
    '''
    if debug:
        t1 = time.time()    
    # Get initial model        
    bins, binx, sol = init_fibermodel(fsize=fsize, bins=bins)

    # Normalize the profiles to 1.
    sol /= np.sum(sol[:-1] * np.diff(binx))
    
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
    ymean = np.mean(Fibers[fib].trace[xlow:xhigh])
    # Define cutouts and fibers
    fibers = Fibers[lowfib:highfib+1]
    x = xgrid[ylow:yhigh,xlow:xhigh].ravel()
    y = ygrid[ylow:yhigh,xlow:xhigh].ravel()
    z = image[ylow:yhigh,xlow:xhigh].ravel()
    # Dummy error
    zerr = np.sqrt(image[ylow:yhigh,xlow:xhigh].ravel())
    
    # selection within cutout
    ycutl=np.repeat(Fibers[lowfib].trace[xlow:xhigh],yhigh-ylow)
    ycuth=np.repeat(Fibers[highfib].trace[xlow:xhigh],yhigh-ylow)
    sel = np.where((y>=ycutl) * (y<=ycuth))[0]    
    # Create empty arrays for the fibermodel weights in each given pixel from 
    Fl = np.zeros((len(y), bins, len(fibers)))
    Pl = np.zeros((len(y),len(fibers)))
    fun = np.zeros((bins,))
    plaw_coeff = np.array([0.0004,0.5,0.15,1.0])
    def plaw(xp, plaw_coeff):
        return plaw_coeff[0] / (plaw_coeff[1] + plaw_coeff[2]
                  * np.power(abs(xp), plaw_coeff[3]))
    i = 0
    for fiber in fibers:
        ytrace = np.repeat([fiber.trace[xlow:xhigh]],(yhigh-ylow))
        ix = y-ytrace
        for j in xrange(bins):
            fun[j] = 1.0
            Fl[:,j,i] = np.interp(ix,binx,fun,left=0.0,right=0.0)
            fun[j] = 0.0
        Pl[:,i] = plaw(ix / 2.5, plaw_coeff)
        i+=1
        
    F = Fl.sum(axis=2)
    P = Pl.sum(axis=1)[sel]    
    
    # Solve for fiber model iteratively with the normalization
    for i in xrange(niter):
        init_model = np.zeros((len(x),len(fibers)))
        flat = np.zeros((len(x),))
        flat_err = np.zeros((len(x),))
        norm = np.zeros((len(fibers),xhigh-xlow))
        normfits = np.zeros((len(x),len(fibers)))
        for j,v in enumerate(np.arange(xlow,xhigh)):
            xsel = np.where(x==v)[0]
            k = 0
            for fiber in fibers:
                init_model[xsel,k] = np.dot(Fl[xsel,:,k],sol) + Pl[xsel,k]
                k+=1
            # Solve for the normalization of the number of fibers
            norm[:,j] = nnls(init_model[xsel,:],z[xsel])[0]
            normfits[xsel,:] = np.ones((len(xsel),1)) * norm[:,j]
            #full_model[xsel] = np.dot(init_model[xsel,:],norm[:,j])
            flat[xsel] = (z[xsel]/((init_model[xsel,:]*normfits[xsel,:]**2)
                         /(init_model[xsel,:]*normfits[xsel,:]).sum(axis=1)
                                                [:,np.newaxis]).sum(axis=1))
            flat_err[xsel] = flat[xsel] / z[xsel] * zerr[xsel]        
     
        try:
            sol1 = nnls(F[sel,1:-1],flat[sel] - P)[0]
        except ValueError:
            print(flat[sel])
            print("nnls for fiber model solution failed on iteration %i" %i)
            sys.exit(1)
        sol = np.hstack([0.0,sol1,0.0])
        sol /= np.sum(sol[:-1] * np.diff(binx))
    return sol
    
    model = np.dot(F,sol) + Pl.sum(axis=1)
    PV = np.array([x,flat_err]).T
    PV = PV[np.argsort(PV[:,0]),:]
    if debug:
        t2 = time.time()
        print("Solution took: %0.3f s" %(t2-t1))        
    if plot:
        fig = plt.figure(figsize=(8,6))
        fibplot = plt.axes([0.1,0.55,0.8,0.4])
        implot = plt.axes([0.1,0.1,0.28,0.4])
        normplot = plt.axes([.55,.1,0.35,0.4])
        fibplot.scatter(x, flat, c=[1.0,0.31,0.31], edgecolor='none', alpha=0.5,
                    s=8)
        ytrace = np.repeat([Fibers[fib].trace[xlow:xhigh]],(yhigh-ylow))
        for i in xrange(len(Fl[0,0,:])):
            lmodel = np.dot(Fl[:,:,i],sol) + Pl[:,i]
            fibplot.scatter(y-ytrace, lmodel, c=[0.41,0.41,0.91],edgecolor='none',
                    alpha=0.5,s=8)
        fibplot.scatter(y-ytrace, model, c=[0.11,0.11,0.91],edgecolor='none',
                    alpha=0.5,s=5)
        fibplot.fill_between(PV[:,0], PV[:,1], -1*PV[:,1],
                         facecolor=[1.0,0.5,.52],edgecolor=[0.9,0.3,.32])
        fibplot.scatter(y-ytrace, flat - model, c=[0.41,0.41,0.41],edgecolor='none',
                    alpha=0.5,s=5)
        fibplot.axis([ylow-ymean, yhigh-ymean, -0.01, 0.22])
        cmap = plt.get_cmap('RdBu_r')
        
        cax = implot.scatter(x[sel],y[sel]-ytrace[sel], c=flat/model, 
                             cmap=cmap, edgecolor='none', s=5, vmin=0.95, 
                             vmax=1.05,marker='s')
        for i in xrange(lowfib,highfib+1):
            yt = Fibers[i].trace[xlow:xhigh]-Fibers[fib].trace[xlow:xhigh]
            implot.plot(np.arange(xlow,xhigh), yt, c=[0.3,0.8,0.32], alpha=0.4)
        implot.axis([xlow, xhigh, ylow-ymean, yhigh-ymean])
        implot.axis('scaled')
        cbaxes = fig.add_axes([0.4, 0.1, 0.03, 0.4]) 
        cb = plt.colorbar(cax, cax = cbaxes)            
        cmap = matplotlib.cm.get_cmap('viridis')
        mx = highfib+1 - lowfib
        for i in xrange(mx):
            normplot.plot(np.arange(xlow+1,xhigh), norm[i,:], 
                          c=cmap((1.*i)/(mx-1)))
        normplot.set_xlim([xlow+1,xhigh])
        fig.savefig("fit_bins_cam%s_%i_%i_%i.png" 
                    %(Fibers[fib].specid, fib, xlow, xhigh), dpi=150)
        plt.close(fig)    
