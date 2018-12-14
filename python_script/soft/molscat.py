#!/usr/bin/env python
'''
splines to extract continuum
Author : Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

from builtins import zip
from builtins import str
import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import scipy.optimize as optimization
import smooth
import argparse


''' fitting lmabda**-4 for molecular scattering'''


def func(x, a, b):
    return a/(x**4) + b


def fromSimu(file):
    values = tb.readtuple(file)
    data = np.array(values).transpose()
    wght = [float(i) for i in data[0]]
    flux = [float(i) for i in data[1]]
    data = np.array([wght, flux]).transpose()
    airmass = float((((os.path.basename(file)).split('_'))[6]).strip('z'))/10.
    Object = 'ref'
    jd = 0.
    parallactic_angle = 0.
    return data, airmass, Object, jd, parallactic_angle


def fromObs(file, colname):
    dict, values, names = tb.readcat(file)
    #print 'Dict : ' ,dict.items()
    latitude = ' '.join(dict.get('LATITUDE'))
    ha = ' '.join(dict.get('HA'))
    dec = dict.get('DEC')
    parallactic_angle = tb.ParallacticAngle(latitude, ha, dec)
    parallactic_angle = parallactic_angle[0]
    print(parallactic_angle)
    jd = ' '.join(dict.get('OBJECT'))
    Object = dict.get('JD')[0]
    airmass = dict.get('AIRMASS')[0]
    data = np.array([values.field('w'), values.field(str(colname))]).transpose()
    #data = np.array([values.field('w'), values.field('psf_gauss_flux')]).transpose()
    return data, airmass, Object, jd, parallactic_angle


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum "

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="analyze spectrum")
    parser.add_argument('-p', "--plot",
                        help="show control plots",
                        action='store_true')
    parser.add_argument('-f', "--outfile", type=str,
                        help="name of output file",
                        default='rawspectrum.list')
    parser.add_argument('-s', "--seeing", type=str,
                        help="seeing from file",
                        default=None)
    parser.add_argument('-c', "--colname", type=str,
                        help="y-axis name",
                        default='aper_flux')
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
    outfile = args.outfile
    inputype = args.Type #''' Either obs or simu '''
    files = args.files
    seeing = args.seeing
    colname = args.colname

    if ((inputype == 'simu') and (seeing is not None)):
        dict, values, names = tb.readcat(seeing)
        seeing_x = values.field('w')
        seeing_y = values.field('psf_gauss_sigma')
        pix2wght = float(dict.get('PIX2WGTH')[0])

    wmin = 300.
    wmax = 1050.

    outlist = []
    for file in files:
        print('opening ', file)
        if (inputype == 'obs'):
            data, airmass, Object, jd, parallactic_angle = fromObs(file, colname)
        if (inputype == 'simu'):
            data, airmass, Object, jd, parallactic_angle = fromSimu(file)
            # from looking at simu
        data = data[(data[:, 0] >= wmin) & (data[:, 0] <= wmax)]
        knots = [[350, 450]]
        #knots=[[350, 450],[778, 782],[860,880],[1035,1045]]

        ''' 
        Smoothing the atmospheric curve to match the resolution of the observation
        '''
        if ((inputype == 'simu') and (seeing is not None)):
            sm = smooth.smooth(data[:, 0], data[:, 1], seeing_x, seeing_y,
                               pix2wgth=pix2wght,
                               plot=plot)
            data[:, 0], data[:, 1] = sm.smoothCurve()

        continuum_w = []
        continuum_f = []
        for k in knots:
            print('continuum : ', k[0], ' to ', k[1])
            continuum_w.extend(data[:, 0][(data[:, 0] >= k[0]) & (data[:, 0] <= k[1])])
            continuum_f.extend(data[:, 1][(data[:, 0] >= k[0]) & (data[:, 0] <= k[1])])

        if len(continuum_w) <= 0:
            continue

        for i, j in zip(continuum_w, continuum_f):
            print(i, j)

        fit, sfit = optimization.curve_fit(func, continuum_w, continuum_f)

        coef = fit[0]
        cst = fit[1]
        print('fitted coef = ', coef, ' cst = ', cst)

        wght = np.linspace(300, 1050, 400)
        molscat = func(wght, coef, cst)

        if plot:
            pl.figure
            pl.plot(data[:, 0], data[:, 1], 'r^', label='raw flux')
            pl.plot(wght, molscat, 'b', label='fitted molecular scattering')
            pl.plot(continuum_w, continuum_f, 'g+', lw=3, label='fitted interval')
            pl.legend()
            pl.show()

    print(outlist)
    names = ['object', 'jd', 'airmass', 'parallactic', 'ew_w', 'ew_f']
    tb.DumpTuple(names,
                 list(zip(*outlist)),
                 outfile)
    print('writing : ', outfile)
