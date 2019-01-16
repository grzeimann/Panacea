# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 16:11:17 2019

@author: gregz
"""

from tables import open_file
import numpy as np
import pyds9


h5file = open_file('/work/03730/gregz/maverick/test.h5', mode='r')
table = h5file.root.Info.Fibers
spec = np.array(table.cols.spectrum[:])
ds9 = pyds9.DS9()
ds9.set_np2arr(spec)
raw_input('Done Inspecting?')
h5file.close()