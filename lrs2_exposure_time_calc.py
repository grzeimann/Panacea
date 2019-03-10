#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 09:15:13 2019

@author: gregz
"""

import numpy as np
from astropy.io import fits
import pandas as pd

from bokeh.layouts import row, column
from bokeh.models import ColumnDataSource, CustomJS, PrintfTickFormatter
from bokeh.models import HoverTool
from bokeh.models.widgets import Slider, Select
from bokeh.plotting import figure, save, output_file

channel = 'UV'
seeing_def = 1.8
objmag_def = 19.0
sky_kind = 'dark'
exptime_def = 1200.0
area = 45e4
fibarea = 3. / 4. * np.sqrt(3.) * 0.59**2
fibarea = 0.452195
colfrac = 0.52
p = fits.open('lrs2_config/positions.fits')
px = p[0].data[0]
py = p[0].data[1]
channels = ['UV', 'Orange', 'Red', 'Farred']
response_dict = {}
for ch in channels:
    F = fits.open('lrs2_config/response_%s.fits' % ch)
    response_dict[ch+'_wave'] = F[0].data[0] * 1.
    response_dict[ch+'_resp'] = F[0].data[1] * 1.

for ch in channels:
    for sname in ['bright', 'grey', 'dark']:
        F = fits.open('lrs2_config/%s_%s_sky.fits' % (sname, ch))
        response_dict[sname + '_' + ch + '_sky'] = F[0].data[1] * 1.    

def get_sn(chan, _seeing, _objmag, skykind, _exptime):
    n1 = 3.* np.sqrt(2)
    wave_arr = response_dict[chan+'_wave']
    wave = np.mean(wave_arr)
    resp_arr = response_dict[chan+'_resp']
    dw = np.mean(np.diff(wave_arr))
    _objflux = 10**(-0.4 * (_objmag - 23.9)) * 1e-29 * 3e18 / wave**2
    _objflux = _objflux / resp_arr
    cnts = _objflux * _exptime * dw * area * colfrac
    _skyflux = response_dict[skykind + '_' + chan + '_sky'] / resp_arr
    skycnts = _skyflux * _exptime * dw * area * fibarea * colfrac
    d = np.sqrt(px**2 + py**2)
    weights = np.exp(-0.5 * d**2 / (_seeing/2.35)**2)
    weights /= weights.sum()
    mask = np.ones(weights.shape, dtype=bool)
    mask[d > (_seeing/2.35*2.5)] = False
    signal = np.sum(cnts * weights[:, np.newaxis]**2 * mask[:, np.newaxis], axis=0) 
    noise = np.sqrt(n1**2 + cnts*weights[:, np.newaxis]*mask[:, np.newaxis] + skycnts)
    noise = np.sqrt(np.sum(noise**2 * weights[:, np.newaxis] * mask[:, np.newaxis], axis=0))
    return wave_arr, signal / noise, noise

# Set up data
x, y, noise = get_sn(channel, seeing_def, objmag_def, sky_kind, exptime_def)
response_dict['x'] = x
response_dict['y'] = y * 0.
DF = pd.DataFrame = response_dict
source = ColumnDataSource(data=DF)


# Set up plot
plot = figure(plot_height=600, plot_width=800, title="Exposure Time Calculator",
              tools="crosshair, pan, reset, save, wheel_zoom",
              x_range=[x.min(), x.max()], y_range=[0.01, 500],
              y_axis_type="log")

plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6,
          line_color='salmon')
plot.xaxis.major_label_text_font_size = "16pt"
plot.yaxis.major_label_text_font_size = "16pt"
plot.xaxis.axis_label = 'Wavelength'
plot.xaxis.axis_label_text_font_size = "20pt"
plot.yaxis.axis_label = 'S/N per pixel'
plot.yaxis.axis_label_text_font_size = "20pt"
plot.xaxis.major_tick_line_color = "firebrick"
plot.xaxis.major_tick_line_width = 3
plot.xaxis.minor_tick_line_color = "orange"
plot.yaxis.major_tick_line_color = "firebrick"
plot.yaxis.major_tick_line_width = 3
plot.yaxis.minor_tick_line_color = "orange"
plot.yaxis[0].formatter = PrintfTickFormatter(format="%3.2f")

plot.add_tools(HoverTool(
    tooltips=[
        ( 'Wavelength',   '@x{0.2f}'            ),
        ( 'S/N per Pixel',  '@y{0.2f}' ), # use @{ } for field names with spaces
    ],

    formatters={
        'Wavelength'      : 'numeral', # use 'datetime' formatter for 'date' field
        'S/N per Pixel' : 'numeral',   # use 'printf' formatter for 'adj close' field
                                  # use default 'numeral' formatter for other fields
    },

    # display a tooltip whenever the cursor is vertically in line with a glyph
    mode='mouse'
))


code="""
    var data = source.data;
    var chw = channel.value + '_wave';
    var chr = channel.value + '_resp';
    var skystr = sky.value + '_' + channel.value + '_sky';
    var w = data[skystr]
    var a = seeing.value / 2.35;
    var b = objmag.value;
    var posx = px;
    var posy = py;
    var d = 0.0
    iweights = new Array(posx.length).fill(0.0);
    weights = new Array(posx.length).fill(0.0);
    mask = new Array(posx.length).fill(0.0);
    var sum = 0.0
    for (var i = 0; i < posx.length; i++){
        d = Math.sqrt(Math.pow(posx[i], 2) + Math.pow(posy[i], 2));
        iweights[i] = Math.pow(2.71828, -0.5 * Math.pow(d / a, 2));
        sum += iweights[i];
        if (d < (a * 2.5))
            mask[i] = 1.0;
    }
    for (var i = 0; i < iweights.length; i++){
        weights[i] = iweights[i] / sum;
    }
    var k = exptime.value;
    var x = data['x'];
    var y = data['y'];
    var wave = data[chw];
    var resp = data[chr];
    var start = data[chw][0];
    var end = data[chw][x.length - 1];
    var dw = wave[1] - wave[0];
    var area = 45E4;
    var fibarea = 0.452195;
    var colfrac = 0.52;
    var cnt = 0.0;
    var skycnt = 0.0;
    var noise = 0.0;
    var n1 = 0.0;
    var signal = 0.0;
    for (var i = 0; i < x.length; i++) {
        cnt = Math.pow(10, (-0.4 * (b - 23.9))) * 1E-29 * 3E18 / Math.pow(start/2 + end/2, 2);
        cnt = cnt / resp[i];
        cnt = cnt * k * dw * area * colfrac;
        skycnt = w[i] / resp[i]
        skycnt = skycnt * k * dw * area * fibarea * colfrac;
        noise = 0.0;
        signal = 0.0;
        for (var j = 0; j < weights.length; j++){
            signal = signal + cnt * weights[j] * weights[j] * mask[j];
            n1 = Math.sqrt(Math.pow(4.24, 2) + cnt * weights[j] + skycnt);
            noise = noise + Math.pow(n1 * weights[j] * mask[j], 2);
        }
        noise = Math.sqrt(noise);
        y[i] = signal / noise;
        x[i] = wave[i];
    }
    sum = 0.0
    for (var i = 0; i < y.length; i++) {
        sum += y[i]
    }
    var avg = sum / y.length
    pt.x_range.setv({"start": start, "end": end})
    pt.title.text = 'Average S/N ' + avg;
    source.change.emit();
"""

# , x_range=plot.x_range
callback = CustomJS(args=dict(source=source, pt=plot, px=px, py=py), code=code)
# Set up widgets
channel_select = Select(value=channel, title='Channel',
                        options=['UV', 'Orange', 'Red', 'Farred'],
                        callback=callback)
sky_select = Select(value='dark', title='Sky Background',
                        options=['dark', 'grey', 'bright'],
                        callback=callback)
seeing = Slider(title='Seeing (")', value=1.8, start=0.8, end=4.0, step=0.1,
                callback=callback)
objmag = Slider(title="Magnitude (AB)", value=19.0, start=6.0, end=25.0,
                step=0.1, callback=callback)

exptime = Slider(title="Exposure Time (s)", value=1200.0, start=5.0, end=3600.,
                 step=5., callback=callback)
callback.args["channel"] = channel_select
callback.args["sky"] = sky_select
callback.args["seeing"] = seeing
callback.args["objmag"] = objmag
callback.args["exptime"] = exptime
#callback.args["channel"] = channel_select

# Set up layouts and add to document
inputs = column(channel_select, sky_select, seeing, objmag, exptime)

output_file("lrs2_exptime_calc.html", title="Exposure Time Calculator")
save(row(inputs, plot, width=1200))