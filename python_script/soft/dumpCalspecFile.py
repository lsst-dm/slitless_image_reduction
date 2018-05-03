#!/usr/bin/env python 
'''
dump spectrum from stis fitsimage

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os, sys
import astropy.io.fits as pf
import toolbox as tb

    
if __name__ == "__main__":
 
    Image = sys.argv[1]
    fits   = pf.open(Image)
    target = fits[0].header.get('TARGETID')
    tbdata = fits[1].data                  # the table is in extension 1
    cols   = tbdata.columns                # names of the columns
    tbdata = tbdata[tbdata.field(0)<11000] # Select lambda<1100nm
    wavelength = tbdata.field(0)/10        # Convert in nm
    flux       = tbdata.field(1)
    staterror  = tbdata.field(2)
    syserror   = tbdata.field(3)

    tb.DumpTuple(('wavelength', 'flux', 'staterror', 'syserror'),\
                 [wavelength, flux, staterror, syserror], str(target)+'.list')
