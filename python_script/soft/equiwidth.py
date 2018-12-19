#!/usr/bin/env python
'''
give a raw fitsimage, position of object in it,
-> return a spectrum

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''
from __future__ import print_function

import os
import sys
import re
import subprocess
import telinst as instru
import numpy as np
import pylab as pl
import spectrum
import astropy.io.fits as pf
import logging
import toolbox as tb
import croaks
import scipy.interpolate as interp

if __name__ == "__main__":

    spectra = sys.argv[1:]
    total = []
    air = []
    fig = pl.figure(1)
    for spectrum in spectra:
        [wgth, flux], keys = tb.readlist(spectrum, ['pixel', 'flux'])
        airmass = keys['AIRMASS']
        air.append(airmass)
        total.append(sum(flux))
        print('airmass : ', airmass)
        pl.plot(wgth, flux, color='black')
        #pl.plot(wgth, flux*float(airmass), color='red')
        spl = interp.UnivariateSpline(wgth, flux, s=1000.)
        xs = np.linspace(100, 700, 600)
        pl.plot(xs, spl(xs), 'g', lw=2)
        pl.show()
        pl.legend()
    fig.savefig("resp.pdf")

    fig = pl.figure(2)
    pl.plot(air, total, 'r^')

    pl.show()
    pl.clf()

    sed = tb.SED(SED)
    fig = pl.figure(2)
    pl.plot(sed.wavelength, sed.flux, color='black')
    pl.title('SED')
    pl.legend()
    fig.savefig("sed.pdf")

    data = np.recfromtxt(atmosphere)
    fig = pl.figure(0)
    pl.plot(data[:, 0], data[:, 1], color='black')
    pl.title('atmospheric transmission')
    pl.legend()
    fig.savefig("atmo.pdf")

    #blur the model
    #atmo[:,1] = filt.gaussian_filter1d(atmo[:,1],sigma=5.8/2.355) #5.8/2.355 is the correct sigma for this instrument/night (and for all of our data)
    #interpolate onto our spectral grid
    #atmointerp = interp.griddata(atmo[:,0],atmo[:,1],s.wavelengths)

    atmo = interp.griddata(data[:, 0], data[:, 1], sed.wavelength)
    sed.flux = sed.flux * atmo

    sed.flux = interp.griddata(sed.wavelength, sed.flux, wgth)
    trans = flux/sed.flux

    fig = pl.figure(1)
