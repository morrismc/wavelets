#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:19:45 2019

@author: matthew

This script is written to replicate the function fft_meanspec written by Adam
booth in matlab.
as of 9/16, I'm still trying to figure out how to get this code working. When 
I ran it on the Non-landslide DEM it took FOREVER to evaluate and looks like
there are some issues, for isntance it has significantly lower power in the non landslide terrane

"""

def fft_mean_spec(DEM, w, dx, normalize, plots):

    import math
    import numpy as np
    from Hann2d import Hann2D
    import progress
    import matplotlib.pyplot as plt
    from scipy import signal
    from detrend_tp import detrend_tp
    

    DEM[DEM == -9999.0] = np.nan
    [nrows, ncols] = np.shape(DEM) # find the dimension of the DEM
    center = int(w/2)+1 # center of moving window
    tpow = 2**(math.ceil(math.log(w)/math.log(2))) #how much zero padding
    Pmat = np.zeros([tpow, tpow]) # intialize the outputgrid

    #calculate frequency bin size (frequency goes from zero to nyquist) as
    #defined by 1/(2*dz) in tpow/2

    df = 1/(dx*tpow)
#    create a matrix of radial frequencies
    xc = tpow/2
    xc = int(xc)
    yc = (xc)

    x = np.array(list(range(int(tpow))))
    y = np.array(list(range(int(tpow))))
    #An annoying extra step required because of different
    #python indexing
#    x = np.delete(x,(0),axis = 0)
#    y = np.delete(y,(0),axis = 0)

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
#    m = 24
#    n = 24
    counter = 0
    for m in range(center,(nrows-center+1)):
        progress.progress(m,(nrows-center+1 - center),'Doing long job')
        for n in range(center,(ncols-center+1)):

#            This next step creates problems, when I go between matlab and python
#            previous to this step all goes well, so I'm not entirely sure
#            what the issue is.
            win = DEM[(m-center+1):(m+center-1),\
                      (n-center+1):(n+center-1)] 
                                                  

            if np.sum(np.isnan(win)) == 0:
                counter = counter + 1
                
                win = signal.detrend(win, type = 'linear')
#                win = detrend_tp(win)

#        %(Optional) Normalize so data has unit variance:
                if normalize == 1:
                    win = win/np.std(win)


                #Variance of the detrended patch
                win_var = np.var(win)

                #window with Hann Raised cosine windo
                win,_ = Hann2D(win)

 #################################### FFT 2d ####################################
                win = np.fft.fftshift(np.fft.fft2(win,[tpow, tpow]))

                #calculate the Discrete fourier periodogram Ampl^2
                win = win*np.conj(win)/(tpow**4)
                win = win.real #necessary b/c otherwise pythong doesn't drop the imaginary
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


#divide by the total number of times through the loop to get mean
    Pmat = Pmat/counter
    
    Pvec = np.reshape(Pmat, tpow*tpow,1)
    fvec = np.reshape(fmat,tpow*tpow,1)
    fp = np.column_stack([fvec,Pvec])
    fp = fp[fp[:,0].argsort(),]
    fvec = fp[:,0]
    Pvec = fp[:,1]
    
#    Pvec = np.sort(Pvec)
#    fvec = np.sort(fvec)
#    Pvec = Pvec[::-1] #reverse the order of Pvec
#    fp = np.column_stack([fvec,Pvec])
    
#    fvec = fp[:,0]
#    Pvec = fp[:,1]

    plt.figure(1)
    plt.imshow(np.log(Pmat))
    
    plt.figure(2)
    plt.loglog(fvec,Pvec,'.')


    return (Pmat, Pvec, fvec, fmat)
