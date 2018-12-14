#!/usr/bin/env python
'''
divide spectra by a spectrum at a reference airmass and fit AOD

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
import matplotlib
import matplotlib.cm as cm
import scipy.interpolate as interp
from mpl_toolkits.mplot3d import Axes3D
import argparse
from scipy.optimize import curve_fit
import smooth
import glob

matplotlib.rcParams.update({'font.size': 6})

OBJECT = ['lamlep', 'ksi02cet', 'hr1544', 'hip96536', 'hip3142', 'hip26382',
          'hip113896', 'hr9087', 'hip107120', 'hd14943', 'mucol', 'hip8417']


'''
fitting the aerosol optical depth
tau = beta * lambda **(-alpha)
'''


class FIT(object):
    def __init__(self, delta_z,
                 plot=False, **kwargs):
        self.delta_z = delta_z

    '''
    from eq 10
    https://arxiv.org/pdf/1210.2619.pdf
    '''

    def tau(self, x, *p):
        alpha, beta = p
        return 10**(0.4 * self.delta_z * beta * (x/1000)**(-alpha))

    def aerosol(self, x, y):
        p0 = [1., 0.007]
        coef, var_matrix = curve_fit(self.tau, x, y,
                                     p0=p0)
        return coef, var_matrix


def measureAOD(W, Z, F, O, length, refW, refF, WS, FS, refWS, refFS,
               plot=False, **kwargs):
    out = []
    colors = iter(cm.rainbow(np.linspace(0, 1, length)))
    fig, ax1 = pl.subplots(figsize=(7, 5))
    pl.xlabel('wavelength (nm)')
    pl.ylabel('flux ratio')
    title = str(O[0])
    pl.title(title)

    refFS = interp.griddata(refWS, refFS, refW)
    for w, z, f, ws, fs in zip(W, Z, F, WS, FS):
        f = interp.griddata(w, f, refW)
        fs = interp.griddata(ws, fs, refW)
        obs_ratio = f/refF
        sim_ratio = fs/refFS
        ratio = obs_ratio/sim_ratio
        delta_z = median_param-z[0]

        fitf = FIT(delta_z)

        Range = [450., 750.]
        wghtr, ratioC = trim(refW, ratio, Range)
        coef, var_matrix = fitf.aerosol(wghtr, ratioC)

        print()
        alpha = coef[0]
        salpha = np.sqrt(var_matrix[0][0])
        tau = coef[1]
        stau = np.sqrt(var_matrix[1][1])
        print('alpha = ', alpha, ' +/- ', np.sqrt(var_matrix[0][0]))
        print(' tau  = ', tau, ' +/- ', stau)
        print()
        out.append([median_param, z[0], alpha, salpha, tau, stau])
        fit = fitf.tau(wghtr, *[alpha, tau])

        if plot:
            fig2 = pl.figure(2)
            pl.plot(refW, ratio, 'r^',
                    label='obs_z_' + str(z[0]) + ' / obs_z_ref')
            pl.plot(wghtr, fit, 'b',
                    label='fit')
            pl.xlabel('wavelength (nm)')
            pl.ylabel('obs_ratio / sim_ratio')
            pl.legend(loc='upper left')
            pl.show()

        color = next(colors)
        ax1.plot(wghtr, fit,
                 color=color,
                 label='airmass '+str(z[0])+'/'+str(median_param))
        ax1.legend(loc='upper right')
    pl.show()

    return out


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
    return


def getMedian(files, param):
    list = []
    for file in files:
        [wght, flux], keys = tb.readlist(file, ['w', 'aper_flux'])
        list.append(float(keys[str(param)]))
    median = np.median(list)
    print('median is ', median)
    return median


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


def trim(wght, flux, Range):
    out = np.array([wght, flux]).transpose()
    out = out[(out[:, 0] >= Range[0])]
    out = out[out[:, 0] <= Range[1]]
    return out[:, 0], out[:, 1]


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
    parser.add_argument('-d', "--doublecheck",
                        help="run sim no_aerosol against sim_aerosol ",
                        action='store_true')

    parser.add_argument('-t', "--path", type=str,
                        help="simu path",
                        default=None)
    parser.add_argument('-i', "--files", nargs='+', type=str,
                        help="input files",
                        default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = grabargs()
    plot = args.plot
    inputype = 'obs'
    files = args.files
    outfile = args.outfile
    simu_path = args.path
    doublecheck = args.doublecheck

    param = 'AIRMASS'
    median_param = getMedian(files, param)

    '''Some list to store all spectra and simu'''
    refW = []
    refWS = []
    refF = []
    refFS = []
    W = []
    WS = []
    F = []
    FS = []
    O = []
    Z = []

    length = len(files)
    Range = [350., 950.]

    if doublecheck:
        #files = [simu_path + 'best/RT_LS_pp_tp_sa_rt_z2047_wv100_oz40.OUT', simu_path + 'best/RT_LS_pp_tp_sa_rt_z1379_wv100_oz40.OUT']
        files = simu_path + 'best/RT_LS_pp_tp_sa_rt_z*_wv100_oz40.OUT'
        files = (glob.glob(files))
        print(files)

    for file in files:
        print('With aerosols :')
        if doublecheck:
            wght, flux, airmass, Object, exptime = extractData(file, 'simu')
        else:
            wght, flux, airmass, Object, exptime = extractData(file, inputype)

        '''extract the simu at same airmass'''
        print('without aerosol ;')
        filesimu = simu_path + 'no_aerosol/RT_LS_pp_tp_sa_rt_z'+str(int(airmass*1000))+'_wv100_oz40.OUT'
        print(filesimu)
        if not os.path.exists(filesimu):
            print(filesimu, ' not found -> continue')
            continue
        else:
            wghtS, fluxS, airmassS, ObjectS, exptimeS = extractData(filesimu, 'simu')

        '''
        Selected a clean range
        '''
        if(airmass >= 2.1):
            print('airmass is above 2.1')
            continue

        '''trimming spectra'''
        wght, flux = trim(wght, flux, Range)
        wghtS, fluxS = trim(wghtS, fluxS, Range)

        flux = flux/exptime
        print(Object, 'airmass ', airmass, len(wght), wght[0], wght[-1])

        if(float(airmass) == median_param):
            refW = wght
            refF = flux
            refWS = wghtS
            refFS = fluxS

        airmass = np.ones(len(wght))*airmass
        W.append(wght)
        F.append(flux)
        WS.append(wghtS)
        FS.append(fluxS)
        Z.append(airmass)
        O.append(Object)

    ''' Plotting spectra normalized by median airmass spectrum'''
    plot1(W, Z, F, O, length, refW, refF)
    plot1(WS, Z, FS, O, length, refWS, refFS)
    pl.show()

    outlist = measureAOD(W, Z, F, O, length, refW, refF, WS, FS, refWS, refFS,
                         plot=plot)

    names = ['median_param', 'z[0]', 'alpha', 'salpha', 'tau', 'stau']
    tb.DumpTuple(names,
                 zip(*outlist),
                 outfile)
    print('writing : ', outfile)
