#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 09:36:10 2019

@author: matthew
"""

def Hann2D(win):
    
    import math
    import numpy as np

    [ny, nx] = np.shape(win)
    #define the centroid of the matrix
    a = nx/2
    b = ny/2
    [X, Y] = np.meshgrid(np.array(list((range(nx)))),
                    np.array(list(range(ny))))
    #now I need to convert to angular coordinates
    theta = (X==a) * (math.pi/2) + (X != a) * np.arctan2((Y-b),(X-a))
    
    r = np.sqrt((Y-b)**(2) + (X-a)**2) # radial polar coords
    rprime = np.sqrt((a**2)*(b**2)*(b**2*(np.cos(theta))**2 +
                      a**2*np.sin(theta)**2)**(-1))
    
    hanncoeff = (r < rprime) * (0.5*(1+np.cos(np.pi*r/(rprime))))                 
    H = win* hanncoeff
    Wss = np.sum(hanncoeff**2)
    
    return(H, Wss)