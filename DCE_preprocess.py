# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 10:22:30 2019
DCE Preprocess, gets all of the inputs needed for the optimized DC_Eig_par function

INPUT: 
    DEM, cellsize, w,
    
OUTPUT:
    cy
    cx
    cz
@author: matthewmorriss
"""
import numpy as np
def DCE_preprocess(DEM, cellsize, w):
    
    [nrows, ncols] = np.shape(DEM)

    #initiate an empty array same size as dem
    rms = DEM*np.nan
#    rms = np.float32(rms)

#    #compute the directional cosines
    [fx, fy] = np.gradient(DEM, cellsize, cellsize)


    grad = np.sqrt(fx**2 + fy**2)
    asp = np.arctan2(fy, fx)



    grad=np.pi/2-np.arctan(grad) #normal of steepest slope
    asp[asp<np.pi]=asp[asp<np.pi]+[np.pi/2]
    asp[asp<0]=asp[asp<0]+[2*np.pi]

    #spherical to cartesian conversion
    r = 1
    cy = r * np.cos(grad) * np.sin(asp)
    cx = r * np.cos(grad) * np.cos(asp)
    cz = r * np.sin(grad)
    
    return(cx,cy,cz)