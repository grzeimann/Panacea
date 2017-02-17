"""

Notes from playing with this in ipython. following this idea should work

In [15]: import matplotlib.pyplot as mpl

In [16]: import matplotlib

In [17]: matplotlib.interactive(True)

In [18]: fig = mpl.figure()

In [19]: fig.add_subplot(4,3,1)
Out[19]: <matplotlib.axes._subplots.AxesSubplot at 0x117291358>

In [20]: f1 = aplpy.FITSFigure("CoFeS20170202T061644.8_073_sci.fits",figure=fig,subplot=(4,3,1))
WARNING: No WCS information found in header - using pixel coordinates [aplpy.header]

In [21]: f1.show_grayscale()
INFO: Auto-setting vmin to -4.667e+01 [aplpy.core]
INFO: Auto-setting vmax to  4.702e+02 [aplpy.core]

In [22]: f2 = aplpy.FITSFigure("CoFeS20170202T061644.8_085_sci.fits",figure=fig,subplot=(4,3,2))
WARNING: No WCS information found in header - using pixel coordinates [aplpy.header]

In [23]: f2.show_grayscale()
INFO: Auto-setting vmin to -6.725e+01 [aplpy.core]
INFO: Auto-setting vmax to  6.930e+02 [aplpy.core]
"""
