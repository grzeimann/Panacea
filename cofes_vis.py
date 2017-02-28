import time
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from glob import glob
import argparse
from collections import Counter
cmap = plt.get_cmap('Greys')
import matplotlib.gridspec as gridspec

def cofes_plots(filename_array, outfile_name, vmin=-15, vmax=25):
    """
    filename_array is an array-like object that contains the filenames
    of fits files to plot. The output plot will be the shape of the input array.
    
    outfile is the output file name including the extension, such as out.fits.
    
    vmin and vmax set the stretch for the images. They are in units of counts;
    Pixels with vmin and below counts are black. Pixels with vmax and above counts
    are white. 
    """
    filename_array = np.array(filename_array)
    assert filename_array.ndim < 3, "filename_array has more than two dimensions. I can't plot that!"
    assert filename_array.size > 0, "filename_array has size zero. There's nothing there to plot!"
    if filename_array.ndim == 1:
        rows = 1
        cols = filename_array.shape[0]
    else:
        rows = filename_array.shape[0]
        cols = filename_array.shape[1]
    
    fig = plt.figure(figsize=(12,12))
    gs1 = gridspec.GridSpec(rows, cols)
    gs1.update(wspace=0.025, hspace=0.025) # set the spacing between axes. 
    for i,fitsfile in enumerate(filename_array.flatten()):
        #robust against files not existing
        try:
            data = fits.open(fitsfile)[0].data
            ax = plt.subplot(gs1[i])
            ax.imshow(data,vmin=vmin,vmax=vmax,interpolation='nearest',origin='lower',cmap=cmap)
            ax.set_xticks([])
            ax.set_yticks([])

            
        except IOError:
            print(fitsfile, "not found. Skipping...")
        
    fig.savefig(outfile_name)

    
    
def cofes_4x4_plots(prefix="", outfile_name = 'CoFeS_plots.png', vmin=-15, vmax = 25):
    """
    dir is a string containing the directory with the CoFeS files you wish
    to plot. If its the local directory you can leave it as an empty string
    
    outfile_name is a string with the output file name. This will be placed
    in the dir directory
    
    the ifu order is
    073 083 093 103
    074 084 094 104
    075 085 095 105
    076 086 096 106

    """
    t1=time.time()
    ifunums = np.array([['073', '083', '093', '103'],
                        ['074', '084', '094', '104'],
                        ['075', '085', '095', '105'],
                        ['076', '086', '096', '106']])
    cofes_files = glob(prefix + "*.fits")
    assert len(cofes_files) > 0, "There are no fits files that begin with " + "'" + \
                                prefix + "'" + ". Please fix your typo."
    prefixes =  [i.split('_')[0] for i in cofes_files]
    prefixes_dict = Counter(prefixes)
    assert len(prefixes) == prefixes.count(prefixes[0]), "The prefix you specified, " + \
                                "'" + prefix + "'" + ", is not unique. It matches to " + \
                                "".join(i + " and " for i in prefixes_dict.keys()) + \
                                ". Please enter a unique prefix."
    prefix = cofes_files[0].split('_')[0]
    filename_list = []
    for i in ifunums.flatten():
        filename_list.append(prefix + '_' + i + '_sci.fits')
    filename_array = np.array(filename_list)
    filename_array = filename_array.reshape(ifunums.shape[0], ifunums.shape[1])
    cofes_plots(filename_array, outfile_name, vmin, vmax)
    
def main():
    """
    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('prefix', type=str, help='Prefix for the CoFeS files you wish to plot. For example CoFeS20170202T061644.8')
    parser.add_argument('--output', default=None, help='Output file name. By default this is the given prefix.png')
    parser.add_argument('--vmin', type=float, default=-15, help='Sets the lower level of stretch for the output images.')
    parser.add_argument('--vmax', type=float, default=25, help='Sets the upper level of stretch for the output images.')
    args = parser.parse_args()

    if args.output is None:
        args.output = args.prefix + '.png'
    #print(args.prefix, args.output, args.vmin, args.vmax)
    cofes_4x4_plots(prefix = args.prefix, outfile_name = args.output, vmin = args.vmin, vmax = args.vmax)
    
    
if __name__ == '__main__':
    main()    
 