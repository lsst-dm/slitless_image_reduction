#!/usr/bin/env python
'''
divide spectra by a spectrum at a reference airmass

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

from builtins import zip
from builtins import next
from builtins import str
import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import croaks
import matplotlib.cm as cm
import scipy.interpolate as interp
from mpl_toolkits.mplot3d import Axes3D
import argparse

OBJECT = ['lamlep', 'ksi02cet', 'hr1544', 'hip96536', 'hip3142', 'hip26382',
          'hip113896', 'hr9087', 'hip107120', 'hd14943', 'mucol', 'hip8417']


def plot1(W, Z, F, O, length, refW, refF):
    colors = iter(cm.rainbow(np.linspace(0, 1, length)))
    fig, ax1 = pl.subplots(figsize=(7, 5))
    pl.xlabel('wavelength (nm)')
    pl.ylabel('flux ratio')
    title = str(O[0])
    pl.title(title)
    for w, z, f in zip(W, Z, F):
        flux_interp = interp.griddata(w, f, refW)
        color = next(colors)
        ax1.plot(refW, flux_interp/refF,
                 color=color,
                 label='airmass '+str(z[0])+'/'+str(median_param))
        ax1.legend(loc='upper right')
    pl.show()
    return


def getMedian(files, param):
    list = []
    for file in files:
        [wght, flux], keys = tb.readlist(file, ['w', 'aper_flux'])
        list.append(float(keys[str(param)]))
    median = np.median(list)
    for i in list:
        print(i)
    return median


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum and extract EW"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="analyze spectrum")
    parser.add_argument('-p', "--plot",
                        help="show control plots",
                        action='store_true')
    parser.add_argument('-f', "--outfile", type=str,
                        help="name of output file",
                        default='default.list')

    parser.add_argument('-t', "--Type", type=str,
                        help="either obs or simu",
                        default=None)
    parser.add_argument('-i', "--files", nargs='+', type=str,
                        help="input files",
                        default=None)
    args = parser.parse_args()
    return args


def extractData(file, inputype):
    print('reading ', file)
    if (inputype == 'obs'):
        [wght, flux], keys = tb.readlist(file, ['w', 'aper_flux'])
        airmass = keys['AIRMASS']
        exptime = keys['EXPTIME']
        Object = keys['OBJECT']

    if (inputype == 'simu'):
        wght, flux, airmass, Object = fromSimu(file)
        exptime = 1.
    return wght, flux, airmass, Object, exptime


def fromSimu(file):
    values = tb.readtuple(file)
    data = np.array(values).transpose()
    wght = [float(i) for i in data[0]]
    flux = [float(i) for i in data[1]]
    try:
        airmass = float((((os.path.basename(file)).split('_'))[6]).strip('z'))/1000.
    except:
        airmass = 0
    Object = 'ref'
    # temporary fix to multiply libradtran output by cos(sza)
    flux = [x * airmass for x in flux]
    return wght, flux, airmass, Object


if __name__ == "__main__":
    args = grabargs()
    plot = args.plot
    inputype = args.Type #''' Either obs or simu '''
    files = args.files
    outfile = args.outfile

    if (inputype == 'obs'):
        param = 'AIRMASS'
        median_param = getMedian(files, param)

    if (inputype == 'simu'):
        #median_param = 1.500
        median_param = 1.379

    length = len(files)
    colors = iter(cm.rainbow(np.linspace(0, 1, length)))
    Range = [350., 950.]
    refW = []
    refF = []
    W = []
    Z = []
    F = []
    O = []
    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')

    for file in files:
        wght, flux, airmass, Object, exptime = extractData(file, inputype)

        '''
        Selected a clean range
        '''
        if(airmass >= 2.1):
            print('airmass is above 2.1')
            continue

        out = np.array([wght, flux]).transpose()
        out = out[(out[:, 0] >= Range[0])]
        out = out[out[:, 0] <= Range[1]]
        if len(out) == 0:
            print('trouble : ', file)
            continue
        wght = out[:, 0]
        flux = out[:, 1]/exptime
        print(Object, 'airmass ', airmass, len(wght), wght[0], wght[-1])

        if(float(airmass) == median_param):
            refW = wght
            refF = flux

        airmass = np.ones(len(wght))*airmass
        W.append(wght)
        Z.append(airmass)
        F.append(flux)
        O.append(Object)
        ax.scatter(wght, airmass, flux, linestyle='-', linewidths=2)

    pl.xlabel('wavelength (nm)')
    pl.ylabel('airmass')
    fig.savefig("all.pdf")

    ''' Plotting spectra normalized by median airmass spectrum'''
    plot1(W, Z, F, O, length, refW, refF)
