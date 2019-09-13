#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:44:04 2019

This is Code built for one purpose: well a couple purposes
1) to load the landslide DEM
2) load the non landslide DEm
3) call fft mean spec so I can begin debugging that code.
@author: matthew
"""


import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from osgeo import gdal

from matplotlib.cbook import get_sample_data
from matplotlib.colors import LightSource

import rasterio as rio
from rasterio.plot import show
from rasterio.mask import mask

#%% load DEMS

os.chdir(r'U:\GHP\Projects\NSF - Morris Landslides\Code\Developmemnt')

matplotlib.rcParams['figure.figsize'] = (8, 5.5)

dem_path = os.path.join(os.getcwd(),'landslide_dem.txt')
demFLD = gdal.Open(dem_path)
demFLD = np.array(demFLD.GetRasterBand(1).ReadAsArray())

dem_path = os.path.join(os.getcwd(),'dem_no_ls.txt')
demUNfld = gdal.Open(dem_path)
demUNfld = np.array(demUNfld.GetRasterBand(1).ReadAsArray()) #this removes
#geogrpahic references

w = 47 # width of the window
dx = 2.0 # cell size

#
#with rio.open(dem_path,'r+') as ds:
#    arr = ds.read() #read all raster values
#    ds.write(arr)
#    
#dem = rio.open(dem_path)
#dem.bounds
#fix,ax = plt.subplots(figsize = (8,3))
#show(dem,
#     title = "lidar",
#     ax = ax)
#ax.set_axis_off()

#%% call fft_mean_spec

from fft_mean_spec import fft_mean_spec
normalize = 0
plots = 1

fft_mean_spec(demFLD, 47, 2, normalize,plots)