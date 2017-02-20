import matplotlib.pyplot as plt
import aplpy
import numpy as np

def cofes_plots(filename_array, outfile_name):
    """
    filename_array is an array-like object that contains the filenames
    of fits files to plot. The output plot will be the shape of the input array.
    outfile is the output file name including the extension, such as out.fits.
    """
    filename_array = np.array(filename_array)
    assert filename_array.ndim < 3, "filename_array has more than two dimensions. I can't plot that!"
    assert filename_array.size > 0, "filename_array has size zero. There's nothing there to plot!"
    if filename_array.ndim == 1:
        rows = 1
        cols = filename_array.shape[0]
    else:
        rows = filename_array.shape[0]
        cols = filename_array.shape[0]
    
    fig = plt.figure()
    for i,fitsfile in enumerate(filename_array.flatten()):
        #robust against files not existing
        try:
            fitsplot = aplpy.FITSFigure(fitsfile,figure=fig,subplot=(rows,cols,i+1))
            fitsplot.show_grayscale()
        except IOError:
            print(fitsfile, "not found. Skipping...")
    
    fig.savefig(outfile_name)
    
    
    
    
    