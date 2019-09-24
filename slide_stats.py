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

def slide_stats(mapped, identified, bins, low, high)

