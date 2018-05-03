#!/usr/bin/env python 
'''
return plots of spectra

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os, sys, re
import numpy as np
import pylab as pl
import toolbox as tb
import croaks
import scipy.interpolate as interp


if __name__ == "__main__":

    wmin    = float(sys.argv[1])
    wmax    = float(sys.argv[2])
    file    = sys.argv[3]
    outfile = sys.argv[4]

    value, keys = tb.readlist(file, ('pixel', 'w','aper_flux','psf_gauss_flux','psf_gauss_sigma', 'psf_gauss_mu','psf_voigt_flux','psf_voigt_fwhmL', 'psf_voigt_fwhmG'))
    

    second = np.array([value[1], value[3]]).transpose()
    second = second[(second[:,0]>= 2*wmin) & (second[:,0]<= 2*wmax)]

    first = np.array([value[1], value[3]]).transpose()
    first = first[(first[:,0]>= wmin) & (first[:,0]<= wmax)]

    interp = interp.griddata(2*first[:,0], first[:,1], second[:,0])
     
    pl.yscale('log')
    pl.plot(2*first[:,0], first[:,1], label = 'first order')
    pl.plot(second[:,0], second[:,1], label = 'contamination')
    pl.plot(second[:,0], second[:,1]/interp, label = 'fraction')
    pl.legend()
    tb.DumpTuple(('wavelength', 'contamination'),(second[:,0], second[:,1]/interp)
                     , outfile )
    
    pl.show()
