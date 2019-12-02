# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:20:28 2019
The actual work of creating the wavelet and conducting the 2d convolution

@author: matthewmorriss
"""

def conv2_mexh(patch,a,dx):
    import numpy as np
    from scipy import signal
    #Generated the Mexican hat wavelet kerne at wavelet scale a.  The kernal much be large enough for the wavelet to decay to ~0 at the edges. The Mexican hat is proportional to the second derivative of a gaussian
    [X,Y] = np.meshgrid(np.arange(-8*a,8*a),np.arange(-8*a,8*a))
    psi = (1/a)*(2 - (X/a)**2 - (Y/a)**2)*(np.exp(-((X/a)**2 + (Y/a)**2)/2))
    # TO PLOT PSI UNCOMMENT
#    from matplotlib import cm
#    ax = plt.axes(projection = "3d")
#    ax.plot_surface(X, Y, psi,cmap=cm.coolwarm)
#   C = 
    # TRYING TO FIGURE OUT THE MOST EFFICIENT 2D CONVOLUTION
#    start = time.time()    
#    C = (dx**2)*signal.convolve2d(patch, psi,'same')
#    end = time.time()
#    print(end-start)
    
    
#    
#    print()
#    print()
#    start = time.time()
    C = (dx**2)*signal.fftconvolve(patch,psi,'same')
#    end = time.time()
#    print(end-start)
#    start = time.time()  
#    C1 = scipy.ndimage.convolve(patch, psi, mode = 'constant', cval = 0)
#    end = time.time()
#    print(end-start)
    # ADAM BOOTH'S SOLUTION, DIDN'T WORK
#    C2 = np.real(np.fft.ifft2(np.fft.fft2(patch) * np.fft.fft2(psi, s = patch.shape)))
#    C = (dx**2)*signal.fftconvolve(patch,psi,'same')
    #convolve patch with the psi using the 2d convolution, multiplying by dx^2 to approximate the double integral. "same" crops C to the same size as the patch
    
    return(C)