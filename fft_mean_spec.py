#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:19:45 2019

@author: matthew

This script is written to replicate the function fft_meanspec written by Adam
booth in matlab.
"""

def fft_mean_spec(DEM, w, dx, normalize, plots):

    import math
    import numpy as np
    import scipy as scipy
    from scipy import linalg
    from scipy import signal
    import sys
    from Hann2d import Hann2D


    DEM[DEM == -9999.0] = np.nan
    [nrows, ncols] = np.shape(DEM) # find the dimension of the DEM
    center = int(w/2)+1 # center of moving window
    tpow = 2**(math.ceil(math.log(w)/math.log(2))) #how much zero padding
    tpow2 = (int(tpow), int(tpow))
    Pmat = np.zeros(tpow2) # intialize the outputgrid

    #calculate frequency bin size (frequency goes from zero to nyquist) as
    #defined by 1/(2*dz) in tpow/2

    df = 1/(dx*tpow)
#    create a matrix of radial frequencies
    xc = tpow/2 + 1
    yc = xc

    x = np.array(list(range(int(tpow+1))))
    y = np.array(list(range(int(tpow+1))))
    #An annoying extra step required because of different
    #python indexing
    x = np.delete(x,(0),axis = 0)
    y = np.delete(y,(0),axis = 0)

    [cols, rows] = np.meshgrid(x, y)
#     %column and row indices
#%Frequency matrix.  Note that since fmat contains an even number of rows
#%and columns, xc and yc are rounded up in row and column indices (shifted
#%down and to the right).  The first row and first column therefore contain
#%the Nyquist frequency in the x- and y-directions, while the last row and
#%column stop at one bin (df) below the Nyquist:

    fmat = ((df*(rows-yc))**2 +\
                     (df*(cols-xc))**2)**(0.5) #raise to exponent.

    #Do a #2D fft in a moving window of size w x w, summing on the go
#    bar = Bar("Processing outer FFT loop", max = (nrows-center+1))
    # for the test
    m = 24
    n = 24
    counter = 0
    for m in range(center,(nrows-center+1)):
#        bar.next()
        for n in range(center,(ncols-center+1)):

            win = DEM[(m-center+1):(m+center-1),\  # # This step creates problems
                      (n-center+1):(n+center-1)]  # # When I go from matlab to python
                                                  #  # I no longer get the exact same #'s

            if np.sum(np.isnan(win)) == 0:
                counter = counter + 1
                win = signal.detrend(win, type = 'linear')

#        %(Optional) Normalize so data has unit variance:
                if normalize == 1:
                    win = win/np.std(np.asarray(win).reshape(-1))




                win_var = np.var(win)

                #window with Hann Raised cosine windo
 #################################### Hanning Window ####################################

                #%Hanning window step
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
                win = H
#                [win] = Hann2D(win)
            #% Do a 2D fft padding with zeros
            # July 14, 2019 reasonably sure everything before this point works
 #################################### FFT 2d ####################################
                win = np.fft.fftshift(np.fft.fft2(win,int(tpow)))

                #calculate the Discrete fourier periodogram Ampl^2
                win = win*np.conj(win)/(tpow**4)

            #%Set power to zero at the zero frequency (DC).  After windowing
            #%the data, its mean may no longer be zero, so this ensures that
            #%the first-order trend is removed:
                win[xc,yc]= 0
#            %Normalize so that the sum of the periodogram equals the
#            %variance of the detrended local patch.  This corrects for the
#            %reduction in variance caused by the windowing function:
                win = win_var*win/np.sum(win)
#                Sum up Pout each time through loop for averaging later:
                Pmat = Pmat + win

#%Generate sorted freqency and power vectors.  Note: these vectors
#%are redundant and could be reduced in size by half, but as coded below
#%they sum to the variance of the original data:


            bar.finish()

    Pvec = np.reshape(Pmat, tpow*tpow,1)
    fvec = np.reshape(fmat,tpow*tpow,1)
    FP = [fvec, Pvec]
    FP = np.sort(FP, axis = -1)
    fvec = FP[:,1]
    Pvec = FP[:,2]

    return (Pmat, Pvec, fvec, fmat)
