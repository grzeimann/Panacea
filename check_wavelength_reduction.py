#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 08:43:49 2019

@author: gregz
"""

import glob
import numpy as np
import sys

from astropy.io import fits
from astropy.table import Table
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, RangeTool
from bokeh.plotting import figure, save, output_file


palette = ["#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
           "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"]

dates = Table.read(sys.argv[1], format='ascii.no_header')
dates = np.sort(dates['col1'])
side = sys.argv[2]
side_dict = {'uv': [4100., 4150.], 'orange':[5650., 5850.], 
             'red': [7450., 7650.], 'farred': [9100., 9300.]}


p = figure(plot_height=300, plot_width=800,
           toolbar_location=None, x_axis_location="above",
           background_fill_color="#efefef", x_range=(side_dict[side][0],
           side_dict[side][1]),
           y_axis_type="log", y_range=(1., 10**5))

select = figure(title=("Drag the middle and edges of the selection "
                       "box to change the range above"),
                plot_height=130, plot_width=800, y_range=p.y_range,
                y_axis_type="log",
                tools="", toolbar_location=None,
                background_fill_color="#efefef")


fn = []
cnt = 0
for date in dates:
    fns = glob.glob('/work/03946/hetdex/maverick/LRS2/CALS/cal_%s_%s.fits' %
                    (date, side))
    for f in fns:
        cnt1 = cnt % len(palette)
        fn.append([f, palette[cnt1], date])
        cnt += 1

source = []
for f in fn:
    F = fits.open(f[0])
    try:
        wavelength = F['response'].data[0]
        counts = np.median(F['arcspec'].data, axis=0)
    
        source.append(ColumnDataSource(data=dict(wavelength=wavelength, 
                                                 counts=counts)))
        p.line('wavelength', 'counts', source=source[-1], line_color=f[1],
               legend='%s' % f[2])
        p.yaxis.axis_label = 'Counts'
        select.line('wavelength', 'counts', source=source[-1], line_color=f[1])
        select.ygrid.grid_line_color = None
    except:
        print('Could not plot %s' % f)

range_tool = RangeTool(x_range=p.x_range)
range_tool.overlay.fill_color = "navy"
range_tool.overlay.fill_alpha = 0.2
select.add_tools(range_tool)
select.toolbar.active_multi = range_tool
p.legend.location = "top_left"

output_file("image.html", title="image.py example")
save(column(p, select))