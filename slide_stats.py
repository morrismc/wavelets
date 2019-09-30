#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 16:14:53 2019
This is a port to python of the Adam Forte code "slide_stats" which is used to step through a series of cutoff values to classify identified landslides and count the correct and incorrect nodes as compared with the mapped landslides

# # INPUTS # #
• mapped = map of 1's and 0's indicating locations of known landslides
• identified = array of spectral power sums produced with 2D DFT or 2D CWT
• bins = number of cutoff value intervals
• low = lowest cutoff value
• high = highest cutoff value
%
# # RETURNS # #
• k = vector of cutoff values
• cls = correctly identified landslide nodes
• cns = correctly identified non-landslide nodes
• ils = incorrectly identified landslide nodes
• ins = incorrectly identified non-landslide nodes
• numNaNs = number of NaNs in array
@author: Matthew C. Morriss written on 9/24/19
"""
import numpy as np


def slide_stats(mapped, identified, bins, low, high):

    if low == 0:
        print('Min and Max values in identified will be used as low high')

        # find min and max in identified
        low = np.minimum(identified)
        high = np.maximum(identified)

    #size of the array for the identified CWT
    [nrows, ncols] = np.shape(identified)

    mapped[np.isnan(identified)|np.isnan(mapped)] = np.nan
    numNaNs = np.sum(np.isnan(mapped))

    # set step size
    h = (high - low)/bins

    #initialize slidemap
    slidemap = np.zeros(nrows, ncols)
    clss = np.zeros(0,bins)
    cns = np.zeros(0,bins)
    ils = np.zeros(0,bins)
    ins = np.zeros(0,bins)


    #start counting ind in output vectors:
    ind = 0


    for k in np.arange(low,high,h):

        #counts time through the loops
        ind = ind + 1

        #setup all values >= k to 1, all values < k to 0
        slidemap = identified >= k

        clss[ind] = sum(sum(slidemap == 1 & mapped == 1,2),1)
        cns[ind] = sum(sum(slidemap == 0 & mapped == 0,2),1)
        ils[ind] = sum(sum(slidemap > mapped,2),1)
        ins[ind] = sum(sum(slidemap < mapped,2),1)

    k = np.transpose(np.arange(low,high,h))



    return(k, clss, cns, ils, ins, numNaNs)


