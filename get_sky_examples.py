#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:37:24 2019

@author: gregz
"""
import glob

path = '/work/03946/hetdex/maverick/LRS2/[UT,PSU,M]*/multi_201902*uv.fits'
files = glob.glob(path)
print(len(files))