# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 15:43:39 2019
This is a copy of the DC Eigen vector code, however, I am trying to paralelize it to run moore efficiently using Numba and the solution provided through stack overflow.

@author: matthewmorriss
"""
import numpy as np
import numba as nb

#
#@nb.njit()
#def isnan(win):
#    for i in range(win.shape[0]):
#        for j in range(win.shape[1]):
#            if np.isnan(win[i,j]):
#                return True
#    return False



@nb.njit(parallel=True)
def DC_eig_par(DEM,w,cx,cy,cz,eps):
    [nrows, ncols] = np.shape(DEM)
#
#    #initiate an empty array same size as dem
    rms = DEM*np.nan
    rms.astype(np.float32)

    #Compute RMS cycling through the DEM
    nw=(w*2)**2

    for i in nb.prange(w+1,nrows-w):


        for j in range(w+1,(ncols-w)):
#            d1=np.int16(np.linspace(i-w,i+w,11))
#            d2=np.int16(np.linspace(j-w,j+w,11))

            tempx = cx[i-w:i+w,j-w:j+w]
            tx = tempx.flatten()
#            tx=np.reshape(tempx,-1)
            

            tempy = cy[i-w:i+w,j-w:j+w]
            ty = tempy.flatten()
#            ty=np.reshape(tempy,-1)
            
#            tempz = np.empty([10,10], dtype = np.float32)
            tempz = cz[i-w:i+w,j-w:j+w]
            tz = tempz.flatten()
#            tz=np.reshape(tempz,-1)
            
            
            if (np.isnan(np.concatenate((tx,ty,tz)))).sum() == 0:
                T=np.array([[np.sum(tx**2), np.sum(tx*ty), np.sum(tx*tz)],
                             [np.sum(ty*tx), np.sum(ty**2), np.sum(ty*tz)], 
                             [np.sum(tz*tx), np.sum(tz*ty), np.sum(tz**2)]])
                
                
                [Te,_] = np.linalg.eig(T) # this step is a bit different from the matlab version b/c np.eig outputs two values.
                l = (Te/nw)
                l[l<eps] = 0
                rms[i,j] = 1/np.log(l[0]/l[1])
                
                
            else:
                rms[i,j] = np.nan

    return(rms)