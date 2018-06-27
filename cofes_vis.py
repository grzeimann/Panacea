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

def cofes_plots(ifunums, filename_array, outfile_name, vmin=-15, vmax=25):
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

    rows=10
    cols=10

    fig = plt.figure(figsize=(8,8))
#    for i, f in zip(ifunums, filename_array):
    for i in np.arange(10):
        for j in np.arange(1,11):
            ifuname='%02d%01d' % (j,i)
            index = [ind for ind, v in enumerate(ifunums) if ifuname==v]
            if len(index):
                f=filename_array[index[0]]
                ax = plt.subplot(rows, cols, i*cols+j)
                ax.set_xticks([])
                ax.set_yticks([])
                try:
                    data = fits.open(f)[0].data
                    ax.imshow(data,vmin=vmin,vmax=vmax,interpolation='nearest',origin='lower',cmap=cmap)
            
                except IOError:
                    print(f, "not found. Skipping...")
                    ax.imshow(np.zeros((49,49)),vmin=1,vmax=2,interpolation='nearest',origin='lower',cmap=cmap)

    plt.text(-335,-11.0,"01",weight='bold')
    plt.text(-285,-11.0,"02",weight='bold')
    plt.text(-235,-11.0,"03",weight='bold')
    plt.text(-185,-11.0,"04",weight='bold')
    plt.text(-135,-11.0,"05",weight='bold')
    plt.text(-85,-11.0,"06",weight='bold')
    plt.text(-35,-11.0,"07",weight='bold')
    plt.text(15,-11.0,"08",weight='bold')
    plt.text(65,-11.0,"09",weight='bold')
    plt.text(115,-11.0,"10",weight='bold')

    plt.text(-366,470,"0",weight='bold')
    plt.text(-366,420,"1",weight='bold')
    plt.text(-366,370,"2",weight='bold')
    plt.text(-366,320,"3",weight='bold')
    plt.text(-366,270,"4",weight='bold')
    plt.text(-366,220,"5",weight='bold')
    plt.text(-366,170,"6",weight='bold')
    plt.text(-366,120,"7",weight='bold')
    plt.text(-366,70,"8",weight='bold')
    plt.text(-366,20,"9",weight='bold')

    plt.subplots_adjust(wspace=0.025, hspace=0.025)    
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
    ifunums = np.array([['012', '022', '032', '042', '052', '062', '072', '082', '092', '102'],
                        ['013', '023', '033', '043', '053', '063', '073', '083', '093', '103'],
                        ['014', '024', '034', '044', '054', '064', '074', '084', '094', '104'],
                        ['015', '025', '035', '045', '055', '065', '075', '085', '095', '105'],
                        ['016', '026', '036', '046', '056', '066', '076', '086', '096', '106'],
                        ['017', '027', '037', '047', '057', '067', '077', '087', '097', '107'],
                        ['018', '028', '038', '048', '058', '068', '078', '088', '098', '108']])
    ifunums = ['%02d%01d' % (j, i) for i in np.arange(10) for j in np.arange(1, 11)]
    ifunums.remove('010')
    ifunums.remove('011')
    ifunums.remove('012')
    ifunums.remove('017')
    ifunums.remove('018')
    ifunums.remove('019')
    ifunums.remove('020')
    ifunums.remove('029')
    ifunums.remove('090')
    ifunums.remove('099')
    ifunums.remove('100')
    ifunums.remove('101')
    ifunums.remove('102')
    ifunums.remove('107')
    ifunums.remove('108')
    ifunums.remove('109')
    ifunums.remove('054')
    ifunums.remove('055')
    ifunums.remove('056')
    ifunums.remove('064')
    ifunums.remove('065')
    ifunums.remove('066')

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
    for i in ifunums:
        filename_list.append(prefix + '_' + i + '_sci.fits')
    filename_array = np.array(filename_list)
#    filename_array = filename_array.reshape(ifunums.shape[0], ifunums.shape[1])
    cofes_plots(ifunums, filename_array, outfile_name, vmin, vmax)
    
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
 
