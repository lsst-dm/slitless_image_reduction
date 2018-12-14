#!/usr/bin/env python
'''
return plots of regression to zero airmass

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

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


def fromSimu(file):
    values = tb.readtuple(file)
    data = np.array(values).transpose()
    wght = [float(i) for i in data[0]]
    flux = [float(i) for i in data[1]]
    data = np.array([wght, flux]).transpose()
    airmass = float((((os.path.basename(file)).split('_'))[6]).strip('z'))/1000.
    Object = 'ref'
    exptime = 1.
    #because libradtran output is in Irradiance
    flux = [x * airmass for x in flux]
    return wght, flux, airmass, Object, exptime


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


if __name__ == "__main__":
    args = grabargs()
    plot = args.plot
    inputype = args.Type #''' Either obs or simu '''
    files = args.files
    outfile = args.outfile
    colors = iter(cm.rainbow(np.linspace(0, 1, len(files))))
    Range = [350., 950.]
    W = []
    Z = []
    F = []
    O = []
    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')
    i = 1
    for file in files:
        print('reading ', file)
        if (inputype == 'obs'):
            [wght, flux], keys = tb.readlist(file, ['w', 'aper_flux'])
            airmass = keys['AIRMASS']
            exptime = keys['EXPTIME']
            Object = keys['OBJECT']

        if (inputype == 'simu'):
            wght, flux, airmass, Object, exptime = fromSimu(file)

        '''
        Selected a clean range
        '''
        if(airmass >= 3.1):
            print('airmass is above 3.1')
            continue

        out = np.array([wght, flux]).transpose()
        out = out[(out[:, 0] >= 350)]
        out = out[out[:, 0] < 1000]
        if len(out) == 0:
            print('trouble : ', file)
            continue
        wght = out[:, 0]
        flux = out[:, 1]/exptime
        print('exptime : ', exptime)
        print(Object, 'airmass ', i, airmass, len(wght), wght[0], wght[-1])
        airmass = np.ones(len(wght))*airmass

        W.append(wght)
        Z.append(airmass)
        F.append(flux)
        O.append(Object)
        if(i == 1):
            refW = wght[(wght >= Range[0])][(wght <= Range[1])]
        ax.scatter(wght, airmass, flux, linestyle='-', linewidths=2)
        i += 1

    #pl.xlim(300,950)
    pl.xlabel('wavelength (nm)')
    pl.ylabel('airmass')
    #pl.title('Nov 26th 2016, HD14943')

    fig.savefig("regression.pdf")
    #pl.show()
    interpW = []
    interpZ = []
    interpF = []
    for w, z, f in zip(W, Z, F):
        flux_interp = interp.griddata(w, f, refW)
        for i, j in zip(refW, flux_interp):
            print(i, j, np.log(j), np.log10(j))
            interpW.append(i)
            interpZ.append(z[0])
            interpF.append(-2.5*np.log10(j))

    #for w,z,f in zip(interpW,interpZ,interpF):
    #    print w,z,f

    degree = 1
    Wreg = []
    Freg0 = []
    sFreg0 = []
    Freg1 = []
    sFreg1 = []
    for wr in refW:
        wght = []
        airmass = []
        flux = []
        for w, z, f in zip(interpW, interpZ, interpF):
            if(wr == w):
                #print w, z, f
                wght.append(w)
                airmass.append(z)
                flux.append(f)

        clean = np.array([airmass, flux]).transpose()
        clean = (clean[~np.isnan(clean).any(1)])
        airmass = clean[:, 0]
        flux = clean[:, 1]
        fit_coef, cov = np.polyfit(airmass, flux, deg=degree, cov=True)
        sigma1 = np.sqrt(cov[0, 0])
        sigma0 = np.sqrt(cov[-1, -1])
        print(wr, 'slope, intercept = ', fit_coef, ' sigma slope = ', sigma1)
        Wreg.append(wr)
        Freg0.append(fit_coef[-1]) #lowest order last
        sFreg0.append(sigma0)
        Freg1.append(fit_coef[0])
        sFreg1.append(sigma1)
        if plot:
            fig2 = pl.figure()
            points = np.linspace(0, 3.0, 30)
            reg = np.polyval(fit_coef, points)
            pl.plot(airmass, flux, 'r^', label='wght = ' + str(wght[0]))
            pl.plot(points, reg)
            pl.xlabel('airmass')
            pl.ylabel('-2.5*log10(flux) (mag)')
            pl.legend(loc='lower right')
            pl.show()

    fig3 = pl.figure()
    pl.yscale('log')

    wpos = []
    pos = []
    wneg = []
    neg = []
    for i, j in zip(Wreg, Freg0):
        if(j >= 0.):
            wpos.append(i)
            pos.append(j)
        else:
            wneg.append(i)
            neg.append(abs(j))
    pl.plot(wpos, pos, color='red', marker='+', label="pos.")
    pl.plot(wneg, neg, color='blue', marker='+', label="neg.")
    #pl.errorbar(Wreg,Freg0, sFreg0, color='red',  marker='+', linestyle='None', label= "Extrapolation to z=0 ")
    pl.title('0 airmass')
    pl.xlabel('wavelength (nm)')
    pl.ylabel('-2.5*log10(flux) (mag)')
    pl.legend(loc='upper left')
    fig3.savefig("regressionOairmassP0.pdf")

    fig4 = pl.figure()
    pl.errorbar(Wreg, Freg1, sFreg1, color='red', marker='+', linestyle='None', label="regression ")
    pl.xlabel('wavelength (nm)')

    pl.legend(loc='upper right')
    fig4.savefig("regressionOairmassP1.pdf")
    pl.show()

    outlist = [Wreg, Freg0, sFreg0, Freg1, sFreg1]
    names = ['wght', 'intercept', 'sintercept', 'slope', 'sslope']
    tb.DumpTuple(names,
                 outlist,
                 outfile)
    print('writing : ', outfile)
