#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 13:57:05 2019
This is a port of the code written by Matteo Berti

This is the Function that will calculate the Direction  cosine eigen value for a moving window

This technique was first pioneered by McKean and Roering (2004)
@author: matthew
"""


def DC_eigen(DEM, cellsize, w):
    import numpy as np
    import  progress as progress


    DEM[DEM == -32767] = np.nan


    [nrows, ncols] = np.shape(DEM)

    #initiate an empty array same size as dem
    rms = DEM*np.nan

    #compute the directional cosines
    [fx, fy] = np.gradient(DEM, cellsize, cellsize)


    grad = np.sqrt(fx**2 + fy**2)
    asp = np.arctan2(fy, fx)

    del(DEM, fx, fy)


    grad=np.pi/2-np.arctan(grad) #normal of steepest slope
    asp[asp<np.pi]=asp[asp<np.pi]+[np.pi/2]
    asp[asp<0]=asp[asp<0]+[2*np.pi]

    #spherical to cartesian conversion
    r = 1
    cy = r * np.cos(grad) * np.sin(asp)
    cx = r * np.cos(grad) * np.cos(asp)
    cz = r * np.sin(grad)

    del(asp, grad)

    #Compute RMS cycling through the DEM
    nw=(w*2+1)**2
    for i in range(w+1,ncols-w):
        total  = ncols-w
        progress(i,total,'Doing long job')

        for j in range(w+1,(ncols-w)):
            d1=np.int64(np.linspace(i-w,i+w,11))
            d2=np.int64(np.linspace(i-w,i+w,11))

            tx=np.reshape(cx[d1[0]:d1[-1],d2[0]:d2[-1]],-1)
            ty=np.reshape(cy[d1[0]:d1[-1],d2[0]:d2[-1]],-1)
            tz=np.reshape(cz[d1[0]:d1[-1],d2[0]:d2[-1]],-1)

            if np.max(np.isnan(np.concatenate((tx,ty,tz)))) == 0:
                T=np.array([[np.sum(tx**2), np.sum(tx*ty), np.sum(tx*tz)], [np.sum(ty*tx), np.sum(ty**2), np.sum(ty*tz)], [np.sum(tz*tx), np.sum(tz*ty), np.sum(tz**2)]])
                [Te,_] = np.linalg.eig(T) # this step is a bit more complicated because of differnces between matlab and python
                l = np.flip(Te/nw)
                l[l<np.finfo(float).eps] = 0
                rms[i,j] = 1/np.log(l[0]/l[1])
            else:
                rms[i,j] = np.nan
    return(rms)
