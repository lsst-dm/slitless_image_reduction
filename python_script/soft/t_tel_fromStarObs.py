#!/usr/bin/env python
'''
compare regression to zero airmass with SED

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

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
import sed as Sed
import math


OBJECT = ['lamlep', 'ksi02cet', 'hr1544', 'hip96536', 'hip3142', 'hip26382',
          'hip113896', 'hr9087', 'hip107120', 'hd14943', 'mucol', 'hip8417']

'''-------------------------------------'''


class Modeling(object):
    def __init__(self, target,
                 convolve=2.,
                 directory='./',
                 plot=False,
                 pix2wgth=2.1,
                 keys='',
                 seeing_at_wght=[2.3e-06, -3.9e-03, 3.6],
                 **kwargs):
        self.target = target
        self.keys = keys
        self.directory = directory
        self.plot = plot
        self.pix2wgth_coef = pix2wgth
        self.convolve = convolve
        self.seeing_at_wght = seeing_at_wght
        self.Init()

    def Init(self):
        self.sed = Sed.SED(self.target, self.convolve,
                           pix2wgth=self.pix2wgth_coef,
                           plot=self.plot,
                           seeing_at_wght=self.seeing_at_wght)

        self.sed_wavelength, self.sed_flux = self.sed.smoothedSed()


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum and extract EW"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="analyze spectrum")
    parser.add_argument('-p', "--plot",
                        help="show control plots",
                        action='store_true')
    parser.add_argument('-o', "--obs", type=str,
                        help="name of output file",
                        default='rawspectrum.list')
    parser.add_argument('-s', "--star", type=str,
                        help="name of star",
                        default=None)
    parser.add_argument('-i', "--seeing", type=str,
                        help="seeing from file",
                        default='/Users/augustinguyonnet/harvard/prod/lamlep/2017_10_09T07_41_12aper80/rawspectrumm+1.list')
    parser.add_argument('-c', "--calib", type=str,
                        help="calib from CBP",
                        default='/Users/augustinguyonnet/harvard/atmosx_git/atmosx/other/Oct2017QE.list')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = grabargs()
    plot = args.plot
    obs = args.obs
    star = args.star
    calib = args.calib
    seeing = args.seeing

    '''use as reference for seeing VS wght'''
    dict, values, names = tb.readcat(seeing)
    w = values.field('w')
    psf_gauss_sigma = values.field('psf_gauss_sigma')
    pix2wght = float(dict.get('PIX2WGTH')[0])

    print 'Obs regression : ', obs
    [wght, intercept, sintercept], keys = tb.readlist(obs, ['wght', 'intercept', 'sintercept'])

    '''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
    '''to dispersion direction (return poly2 fit parameters)                         '''
    seeing_at_wght = tb.SeeingAtWght(w, psf_gauss_sigma)
    out = [seeing_at_wght[0], seeing_at_wght[1], seeing_at_wght[2]]

    '''Loading SED'''
    mod = Modeling(star,
                   plot=plot,
                   pix2wgth=pix2wght,
                   seeing_at_wght=seeing_at_wght)

    '''Retrieving smoothed SED'''
    sed_wght = mod.sed_wavelength
    sed_flux = mod.sed_flux

    ''' Converting mag back to flux'''
    intercept = [math.pow((x/(-2.5)), 10) for x in intercept]
    sintercept = [math.pow((x/(-2.5)), 10) for x in sintercept]

    norm = max(intercept)
    intercept = np.array(intercept) / norm
    norm = max(sintercept)
    sintercept = np.array(sintercept) / norm
    norm = max(sed_flux)
    sed_flux = np.array(sed_flux) / norm

    '''Looking at SED and obs'''
    fig = pl.figure()
    pl.plot(sed_wght, sed_flux, label='SED')
    pl.errorbar(wght, intercept, sintercept, label='Obs')
    pl.xlabel('wavelength (nm)')
    pl.ylabel('SED (norm at max)')
    pl.legend()

    intercept = interp.griddata(wght, intercept, sed_wght)
    sintercept = interp.griddata(wght, sintercept, sed_wght)

    '''Determination of telescope calibration '''
    t_tel = intercept/sed_flux
    norm = np.nanmax(t_tel)
    t_tel = np.array(t_tel)/norm

    print 'Telescope reference calibration from : ', calib
    [wght2, qe], keys = tb.readlist(calib, ['wght', 'qe'])
    norm = np.nanmax(qe)
    qe = np.array(qe)/norm

    fig = pl.figure()
    pl.plot(wght2, qe, 'r', label='telescope throughput from CBP')
    pl.plot(sed_wght, t_tel, 'b', label='telescope throughput from (linear regression)/SED')
    pl.xlabel('wavelength (nm)')
    pl.ylabel('Throughput (norm at max)')
    pl.legend()
    pl.show()
