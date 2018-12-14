#!/usr/bin/env python
'''
return plots of spectra

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import scipy.interpolate as interp

if __name__ == "__main__":
    if(len(sys.argv) < 2):
        usage = 'give a list of spectrum with pixel and flux columns'
        sys.exit(usage)

    spectra = sys.argv[1:]
    total = []
    air = []
    fig = pl.figure(1)
    for spectrum in spectra:
        [wgth, flux], keys = tb.readlist(spectrum, ['pixel', 'flux'])
        airmass = keys['AIRMASS']
        air.append(airmass)
        total.append(sum(flux))
        print 'airmass : ', airmass
        pl.plot(wgth, flux, color='black')
    fig.savefig("resp.pdf")
    fig = pl.figure(2)
    pl.plot(air, total, 'r^')
    pl.show()
    pl.clf()
