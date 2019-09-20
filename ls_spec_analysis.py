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

import time
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from osgeo import gdal

#%% load DEMS

os.chdir(r'U:\GHP\Projects\NSF - Morris Landslides\Code\Developmemnt')

matplotlib.rcParams['figure.figsize'] = (8, 5.5)

dem_path = os.path.join(os.getcwd(),'ls_patch.tif')
demFLD = gdal.Open(dem_path)
demFLD = np.array(demFLD.GetRasterBand(1).ReadAsArray())

dem_path = os.path.join(os.getcwd(),'no_ls_patch.tif')
demUNfld = gdal.Open(dem_path)
demUNfld = np.array(demUNfld.GetRasterBand(1).ReadAsArray()) #this removes

dem_path = os.path.join(os.getcwd(),'snowBasinFull1.tif')
dem = gdal.Open(dem_path)
dem = np.array(dem.GetRasterBand(1).ReadAsArray()) #this removes
#geogrpahic references
demUNfld[demUNfld == -32767] =np.nan

w = 47 # width of the window, must be odd
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

start = time.time()
[Vdftave_fld, Vdftvec_fld, fvec, freqmat] = fft_mean_spec(demFLD, 47, 0.5, normalize,plots)

[Vdftave_unfld, Vdftvec_unfld, fvec, freqmat] = fft_mean_spec(demUNfld, 47, 0.5, normalize,plots)

from fft_normPower import fft_normPower
#calculate normalized fourier
[Vdft_norm, Vdftvec_norm,fvecN] = fft_normPower(Vdftave_fld, Vdftave_unfld, freqmat, plots)
end = time.time()
print(end-start)
#%% Try plotting a few things
#fig,(ax1,ax2) = plt.subplots(1,2)

plt.loglog(fvec,Vdftvec_fld,'r.', label = 'Landslide Terrane')
plt.loglog(fvec,Vdftvec_unfld,'b.', label = 'Non-landslide')
plt.xlabel('Wavelength')
plt.ylabel('Spectral Power')
plt.title('Spectra for landslide and non-landslide terrane')
plt.legend()
plt.xlim(0.03, 1)
plt.ylim(10e-11, 10e-1)
plt.show()


plt.semilogx(fvecN,Vdftvec_norm,'.')
plt.xlabel('Wavelength')
plt.ylabel('Normalized Power')
plt.title('Normalized spectral power')



#%% now use wavelets
from conv2_mexh_var import conv2_mexh_var
dx = 0.5
scales = np.exp(np.linspace(0,2.2,20))
[Vcwt_fld, frq, wave] = conv2_mexh_var(demFLD, scales, dx)
[Vcwt_unfld, _, _] = conv2_mexh_var(demUNfld, scales, dx)



Vcwt_fld = np.transpose(Vcwt_fld)
Vcwt_unfld = np.transpose(Vcwt_unfld)
Vcwt_norm = Vcwt_fld/Vcwt_unfld
#%% plot results from the wavelet transform
fig1, ax1 = plt.subplots()
ax1.loglog(frq, Vcwt_fld,'s')
ax1.loglog(frq,Vcwt_unfld,'v')


fig2, ax2 = plt.subplots()
ax2.semilogx(frq, Vcwt_norm)

#%% Compute the wavelt coefficients. 
#This is the step where some work is involved in solving for the wavelet scale needed to filter out the wavelenths of interest
# W = 2pi*dx*s/(sqrt(5/2))
# this equation can be rearranged to solve for s which is the input in these next bits
from conv2_mexh2 import conv2_mexh2

[C2,_,_] = conv2_mexh2(dem, 1.5,dx)
[C3,_,_] = conv2_mexh2(dem, 2.0,dx)
[C4,_,_] = conv2_mexh2(dem, 2.5,dx)
[C5,_,_] = conv2_mexh2(dem, 3.0,dx)
[C6,_,_] = conv2_mexh2(dem, 3.5,dx)

#%% Square and sum wavelet coefficients in quadrature

Vcwtsum = C2**2 + C3**2 + C4**2 + C5**2 + C6**2

fig1, ax1= plt.subplots()
plt.imshow(np.log(Vcwtsum))
ax1.title('CWT Spectral Power Summer')

radius = 25
from smooth2 import smooth2
[Vcwt_smooth,_] = smooth2(Vcwtsum,radius)

#fig2, ax2 = plt.imshow()
#ax2.title('CWT Power Sum (smoothed)')
#ax2.imshow(np.log(Vcwt_smooth))
