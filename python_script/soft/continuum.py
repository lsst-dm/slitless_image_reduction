#!/usr/bin/env python
'''
splines to extract continuum
Author : Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

from builtins import zip
import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import scipy.interpolate as interp
import itertools
import smooth
import argparse


'''Extraction of EW from template input'''


def equWidth(wavelength, flux, continuum, list):
    data = np.array([wavelength, flux, continuum]).transpose()
    cut = data[(data[:, 0] >= list[0]) & (data[:, 0] <= list[1])]
    spacing = cut[1, 0] - cut[0, 0]
    equ_w = np.trapz(1 - cut[:, 1] / cut[:, 2], dx=spacing)
    #Sum=0
    #for i,j,k in zip(cut[:,0],cut[:,1], cut[:,2]):
    #    print i,j,k
    #    Sum += 1- j/k
    #    print "ew = " , Sum
    segment = (list[0]+list[1])/2.
    print('pos, EW : ', segment, equ_w)
    return segment, equ_w


'''Extraction of EW, fitting 2-d polynom on edges'''


def equWidth2(wavelength, flux, list, plot=False, **kwargs):
    extent = 10.
    data = np.array([wavelength, flux]).transpose()
    #edge1   = data[(data[:,0]>= list[0]-extent) &  (data[:,0]<= list[0])]
    #edge2   = data[(data[:,0]>= list[1]) &  (data[:,0]<= list[1]+extent)]
    edge1 = data[(data[:, 0] >= list[1][0]) & (data[:, 0] <= list[1][1])] #going inward
    edge2 = data[(data[:, 0] >= list[2][0]) & (data[:, 0] <= list[2][1])] #going inward

    if ((len(edge1) == 0) or (len(edge2) == 0)):
        return 0, 0

    fit1 = np.polyfit(edge1[:, 0], edge1[:, 1], deg=2)
    fit2 = np.polyfit(edge2[:, 0], edge2[:, 1], deg=2)

    e1 = np.polyval(fit1, list[1][1])
    e2 = np.polyval(fit2, list[2][0])

    fit = np.polyfit([list[1][1], list[2][0]], [e1, e2], deg=1)
    wght = np.linspace(list[1][1], list[2][0], int(list[2][0]-list[1][1]))
    continuum = np.polyval(fit, wght)

    signal = interp.griddata(wavelength, flux, wght)

    spacing = wght[1]-wght[0]
    equ_w = np.trapz(1 - signal / continuum, dx=spacing)

    if plot:
        pl.figure()
        pl.plot(wavelength, flux, label='signal')
        pl.plot(edge1[:, 0], edge1[:, 1], 'r^', label='edge1')
        pl.plot(edge2[:, 0], edge2[:, 1], 'r^', label='edge2')

        points1 = np.linspace(list[1][0], list[1][1], 20)
        out1 = np.polyval(fit1, points1)
        pl.plot(points1, out1, linewidth=1, color='k')
        points2 = np.linspace(list[2][0], list[2][1], 20)
        out2 = np.polyval(fit2, points2)
        pl.plot(points2, out2, linewidth=1, color='k')
        pl.plot(wght, continuum, linewidth=1, color='g')
        pl.legend()

        pl.figure()
        pl.plot(wght, (signal) / continuum, label='signal/continuum')

        pl.show()

    segment = (list[2][0]+list[1][1])/2.
    return segment, equ_w


def fromSimu(file):
    values = tb.readtuple(file)
    data = np.array(values).transpose()
    wght = [float(i) for i in data[0]]
    flux = [float(i) for i in data[1]]
    data = np.array([wght, flux]).transpose()
    try:
        airmass = float((((os.path.basename(file)).split('_'))[6]).strip('z'))/1000.
    except:
        airmass = 0

    try:
        pwv = float((((os.path.basename(file)).split('_'))[7]).strip('wv'))/100.
    except:
        pwv = 0
    Object = 'ref'
    jd = 0.
    parallactic_angle = 0.
    return data, airmass, Object, jd, parallactic_angle, pwv


def fromObs(file):
    dict, values, names = tb.readcat(file)
    meanseeing = meanSeeing(values.field('w'), values.field('psf_gauss_sigma'))
    #print 'Dict : ' ,dict.items()
    latitude = ' '.join(dict.get('LATITUDE'))
    ha = ' '.join(dict.get('HA'))
    dec = dict.get('DEC')
    parallactic_angle = tb.ParallacticAngle(latitude, ha, dec)
    parallactic_angle = parallactic_angle[0]
    print('parallactic_angle ', parallactic_angle)
    Object = ' '.join(dict.get('OBJECT'))
    jd = dict.get('JD')[0]
    airmass = dict.get('AIRMASS')[0]
    #data = np.array([values.field('w'), values.field('aper_flux')]).transpose()
    data = np.array([values.field('w'), values.field('psf_gauss_flux')]).transpose()
    print('flux estimator is Gauss PSF')
    return data, airmass, Object, jd, parallactic_angle, meanseeing


#knots=[150, 200, 220, 250, 300, 420, 450, 550, 600, 700, 800, 880, 1000, 1100]
#knots=[680, 730, 795, 850, 880, 1050, 1100]


EW = [['o2', [735, 745], [780, 785]],
      ['pwv', [850, 870], [1030, 1050]]]


def meanSeeing(wght, sigma):
    data = np.array([wght, sigma]).transpose()
    seeing = data[:, 1][(data[:, 0] >= EW[1][1][1]) & (data[:, 0] <= EW[1][2][1])]
    return np.mean(seeing)


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
                        default='regression.list')
    parser.add_argument('-s', "--seeing", type=str,
                        help="seeing from file",
                        default=None)
    parser.add_argument('-l', "--select", type=str,
                        help="select subset",
                        default=None)
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
    select = args.select

    if ((inputype == 'simu') and (seeing is not None)):
        dict, values, names = tb.readcat(seeing)
        seeing_x = values.field('w')
        seeing_y = values.field('psf_gauss_sigma')
        pix2wght = float(dict.get('PIX2WGTH')[0])
        mean_seeing = meanSeeing(seeing_x, seeing_y)
        print('MEAN SEEING ', mean_seeing)

    wmin = 680.
    wmax = 1200.
    width = 2.

    outlist = []
    for file in files:
        print('opening ', file)
        if (inputype == 'obs'):
            data, airmass, Object, jd, parallactic_angle, mean_seeing = fromObs(file)
            water_knots = [[680, 743], [800, 870], [1030, 1050]] # from looking at obs
            order = 1
            pwv = 0.
        if (inputype == 'simu'):
            data, airmass, Object, jd, parallactic_angle, pwv = fromSimu(file)
            water_knots = [[680, 685], [710, 712], [740, 755], [845, 880], [1000, 1030]]
            # from looking at simu
            order = 2
        data = data[(data[:, 0] >= wmin) & (data[:, 0] <= wmax)]

        ''' 
        Smoothing the atmospheric curve to match the resolution of the observation
        '''
        if ((inputype == 'simu') and (seeing is not None)):
            sm = smooth.smooth(data[:, 0], data[:, 1], seeing_x, seeing_y,
                               pix2wgth=pix2wght,
                               plot=plot)
            data[:, 0], data[:, 1] = sm.smoothCurve()
        if seeing is None and inputype is 'simu':
            mean_seeing = 0.
        '''
        Select  subset of the data
        '''
        if (select == 'jd'):
            print('jd, airmass, meanSeeing : ', jd, airmass, mean_seeing)
            #if((float(jd) < 2458037.5) or (float(jd)>2458039.8)):
            if((float(jd) < 2458038.5) or (float(jd) > 2458039.)): #night 2
                print('wrong period for ', file)
                continue

        '''
        Select the continuum
        '''
        continuum_w = []
        continuum_f = []
        for k in water_knots:
            print('continuum : ', k[0], ' to ', k[1])
            #continuum_w.extend(data[:,0][(data[:,0]>= k-width) & (data[:,0]<= k+width)])
            #continuum_f.extend(data[:,1][(data[:,0]>= k-width) & (data[:,0]<= k+width)])
            continuum_w.extend(data[:, 0][(data[:, 0] >= k[0]) & (data[:, 0] <= k[1])])
            continuum_f.extend(data[:, 1][(data[:, 0] >= k[0]) & (data[:, 0] <= k[1])])

        if len(continuum_w) < 5:
            continue

        #for i,j in zip(continuum_w, continuum_f):
        #    print i,j

        spl = interp.UnivariateSpline(continuum_w, continuum_f, k=order, s=0)
        #interp.spline(data[:,0], data[:,1], knots, order=3, kind='smoothest', conds=None)

        wght = np.linspace(680, 1050, 370)
        continuum = spl(wght)
        signal = interp.griddata(data[:, 0], data[:, 1], wght)
        edges = []
        for item in EW:
            edges.append(item[1][1])
            edges.append(item[2][0])
        for i in EW:
            print('EW : ', i)
            #ew_w, ew_f = equWidth(wght, signal, continuum, i) #abrupt version
            ew_w, ew_f = equWidth2(data[:, 0], data[:, 1], i, plot=plot)          #fitted edges
            print(Object, jd, airmass, ew_w, ew_f, mean_seeing)

            outlist.append([Object, jd, airmass, parallactic_angle, ew_w, ew_f, mean_seeing, pwv])

        if plot:
            pl.figure
            pl.plot(wght, continuum, 'g', lw=1, label='continuum')
            pl.plot(data[:, 0], data[:, 1], 'r^', label='raw flux')
            pl.plot(continuum_w, continuum_f, 'bo', label='spline input')
            pl.legend()

            pl.figure()
            pl.plot(edges, np.ones(len(EW)*2), 'r|')
            pl.plot(wght, (signal) / continuum, label='signal/continuum')
            pl.legend()
            pl.show()

    names = ['object', 'jd', 'airmass', 'parallactic', 'ew_w', 'ew_f', 'mean_seeing', 'pwv']
    tb.DumpTuple(names,
                 list(zip(*outlist)),
                 outfile)
    print('writing : ', outfile)
