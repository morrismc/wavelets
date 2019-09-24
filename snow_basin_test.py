# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 13:03:17 2019
This is a script now to test the entire Snowbasin area using the constraints I weas able to develop using the ls_spec_analysis.py

I have been able to constrain the scale of landslides between 13 - 3 scale for 0.5 m data. 
Now I will attempt to load a larger DEM of the area and see how it performs with a wavelet powersum approach, This should highlight areas of mapped landslides.
@author: matthewmorriss
"""

import time
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from osgeo import gdal


#%%
os.chdir(r'U:\GHP\Projects\NSF - Morris Landslides\Code\Developmemnt')

matplotlib.rcParams['figure.figsize'] = (8, 5.5)

dem_path = os.path.join(os.getcwd(),'sb_less_steep.tif')
sbDEM = gdal.Open(dem_path)
sbDEM = np.array(demFLD.GetRasterBand(1).ReadAsArray())


sbDEM[sbDEM == -9999.0] = np.nan
sbDEM[sbDEM == -32767] =np.nan


#View hillshade of the DEM

from hillshade import hillshade
hs_array = hillshade(demFLD)
plt.imshow(hs_array,cmap='Greys')
plt.show()



#%% Now convolve landscape with the MexH at the appropriate scales


from conv2_mexh2 import conv2_mexh2

[C2,_,_] = conv2_mexh2(sbDEM, 13,dx)
[C3,_,_] = conv2_mexh2(sbDEM, 11,dx)
[C4,_,_] = conv2_mexh2(sbDEM, 9,dx)
[C5,_,_] = conv2_mexh2(sbDEM, 7,dx)
[C6,_,_] = conv2_mexh2(sbDEM, 5,dx)
[C7,_,_] = conv2_mexh2(sbDEM, 4,dx)
[C8,_,_] = conv2_mexh2(sbDEM, 3,dx)

#%% Square and sum wavelet coefficients in quadrature

Vcwtsum = C2**2 + C3**2 + C4**2 + C5**2 + C6**2 + C7**2 + C8**2

fig1, ax1= plt.subplots()
plt.imshow(np.log(Vcwtsum))
ax1.set_title('CWT Spectral Power Sum')

radius = 25
from smooth2 import smooth2
[Vcwt_smooth,_] = smooth2(Vcwtsum,radius)

fig2, ax2 = plt.imshow()
ax2.set_title('CWT Power Sum (smoothed)')
ax2.imshow(np.log(Vcwt_smooth))


