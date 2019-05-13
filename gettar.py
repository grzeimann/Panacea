#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 09:06:24 2019

@author: gregz
"""

import os.path as op
import glob
import tarfile
from astropy.table import Table

tarfolder = op.join('/work/03946/hetdex/maverick/201*', 'virus', 
                    "{:s}00000*.tar".format('virus'))
tarfolders = sorted(glob.glob(tarfolder))
filenames = []
kind = []
for tarfolder in tarfolders:
    print(tarfolder)
    T = tarfile.open(tarfolder, 'r')
    flag = True
    while flag:
        a = T.next()
        try:
            name = a.name
        except:
            flag = False
        if name[-5:] == '.fits':
            filenames.append(tarfolder)
            kind.append(name[-8:-5])
            flag = False
t = Table([filenames, kind], names=['Filename', 'Kind'])
t.write('alltar.txt', format='ascii.fixed_width_two_line')