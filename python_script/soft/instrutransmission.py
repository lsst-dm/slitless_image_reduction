#!/usr/bin/env python 
'''
give a raw fitsimage, position of object in it,
-> return a spectrum

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os, sys, re
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
import sed as Sed


def readtxt(inputfile, dic=[]):
    data = croaks.NTuple.fromtxt(inputfile)
  
    value =[]
    for name in dic:
        value.append(np.array(data[:][name]))
    return value, data.keys


if __name__ == "__main__":
    narg = len(sys.argv)
    if narg<2 :
        print "process.py [-i -j fitsimage] of [fitsimages]"
        print "If keword is none, print whole header"
        print
    k = 1
    while( k<narg ):
        if( sys.argv[k] == "-r" ): 
            k += 1
            raw_spectrum = sys.argv[k]
            k += 1
        elif( sys.argv[k] == "-s" ):
            k += 1
            SED     = sys.argv[k]
            k += 1
        elif( sys.argv[k] == "-a" ):
            k += 1
            atmosphere = sys.argv[k]
            k += 1


    sed = Sed.SED(SED)
    fig = pl.figure(2)
    pl.plot(sed.wavelength, sed.flux, color='black')
    pl.title('SED')
    pl.legend()
    fig.savefig("sed.pdf")

    
    data    = np.recfromtxt(atmosphere)
    fig = pl.figure(0)
    pl.plot(data[:,0], data[:,1], color='black')
    pl.title('atmospheric transmission')
    pl.legend()
    fig.savefig("atmo.pdf")


    
    atmo     = interp.griddata(data[:,0], data[:,1], sed.wavelength)
    sed.flux = sed.flux * atmo

    [ wgth,flux], keys  = readtxt(raw_spectrum, ['w', 'rawflux'])
    airmass = keys['AIRMASS']
    print 'airmass : ', airmass

    sed.flux = interp.griddata(sed.wavelength, sed.flux, wgth)
    trans   = flux/sed.flux
    
    fig = pl.figure(1)
    pl.plot(wgth, trans, color='black', label='Instrumental transmissison' )
    pl.legend()
    fig.savefig("resp.pdf")
    pl.show()
    pl.clf()
