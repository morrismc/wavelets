#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:09:53 2019

@author: matthew
"""

#Test  of the DC Eigen

import time
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from osgeo import gdal

os.chdir(r'/Users/matthew/Dropbox/Post_Phd_jobs/Utah/NSF_Internship/Data')
dem_path = os.path.join(os.getcwd(),'ls_patch_4.tif')
DEM = gdal.Open(dem_path)
DEM = np.array(DEM.GetRasterBand(1).ReadAsArray())


os.chdir(r'/Users/matthew/Documents/GitHub/wavelets')

cellsize = 0.5
w = 5


#%% try eigen function

from DC_eigen import DC_eigen

r = DC_eigen(DEM, cellsize,w)

