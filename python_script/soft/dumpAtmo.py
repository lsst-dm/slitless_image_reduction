#!/usr/bin/env python 
'''
dump a atmospheric trnasmission profile

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os, sys
import toolbox as tb
import scipy.ndimage.filters as filt
    
if __name__ == "__main__":
    wavelength, flux = tb.AtmoTransmission()
    flux             = filt.gaussian_filter1d(flux, sigma=2.46)
    tb.DumpTuple(('wavelength', 'flux'),[wavelength, flux], 'atmo.list')
