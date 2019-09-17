# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 14:42:34 2019
Lsplane.py Least square regression plane

INPUT
X   array [X Y Z] where x = vector of x coordinates, y = y coords, z = z coords

output
x0 =centroid of the data = point on the best fit plane
dim 3x1

a = direction cosines of the normal to the best-fit plane
dim 3x1

optional, removed for this code

d = residuals
m x 1

normd norm of the residual errors
dim 1x1

@author: matthewmorriss
"""

def lsplane(mat):
    import numpy as np
    import sys
    import scipy
    from scipy import linalg
    #check for the number of data points
    [_,m] = np.shape(mat)
    
    if m <3:
        print('Error' )
        sys.exit()
        
    #calculate centroid
    x0 = np.transpose(np.mean(mat, axis = 0))
    
    #form matrix A of translated points
    #Important note about this product, there are decimals and fractions due to the issues surrounding indexing in python... not sure if it'll be an issue or not.
    A = np.column_stack([(mat[:,0]-x0[0]),\
                        (mat[:,1]-x0[1]),\
                        (mat[:,2]-x0[2])])
        
    #calculate the SVD of A
    _,S,Vh = scipy.linalg.svd(A, overwrite_a = True, check_finite = False)
    V = Vh.T.conj()
    S = np.diag(S)
    
    #find the smallest singular value in S and extract from V the corresponding right singular vector
    [_, i] = [min(np.diag(S)),np.argmin(np.diag(S))]
    a = V[:,i]
    

    
    return(x0, a,)