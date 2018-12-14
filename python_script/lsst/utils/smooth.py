#!/usr/bin/env python

import os
import sys
import re
import astropy.io.fits as pf
import numpy as np
import pylab as pl
import re
import logging
import scipy.ndimage.filters as filt
import scipy.interpolate as interp
import toolbox as tb

'''
Smooth data array with seeing array
'''

# ===========================================================


class smooth(Exception):
    pass
# ===========================================================


'''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
'''to dispersion direction (return poly2 fit parameters)                         '''


'''Match target name with its calibrated spectrum'''


class smooth(object):
    def __init__(self, wgth, data, seeing_x, seeing_y,
                 unit='nm',
                 pix2wgth=2.,
                 plot=False, **kwargs):
        self.plot = plot
        self.wgth = wgth
        self.data = data
        self.pix2wgth_coef = pix2wgth
        self.profile_param = self.profile(seeing_x, seeing_y)
        self.init(unit)

    def init(self, unit):
        print 'unit is in ', unit
        self.step = float(self.wgth[1] - self.wgth[0])
        print 'step = ', self.step, unit

    '''fit a parabolic profile from sigma(in pixel) VS wght data'''

    def profile(self, seeing_x, seeing_y):
        clean = np.array([seeing_x, seeing_y]).transpose()
        clean = clean[(clean[:, 0] > 700) & (clean[:, 0] < 1000)]
        profile_param = np.polyfit(clean[:, 0], clean[:, 1], deg=2)
        return profile_param

    ''' Geometric distorsion as a function of angle (wavelength)'''
    ''' due to defocus of a non collimated beam                 '''

    def Resolution(self, item, method='local'):
        local_seeing = np.polyval(self.profile_param, item)
        local_seeing = local_seeing * self.pix2wgth_coef / self.step
        # convert in nm, then to the data bin size
        return local_seeing

    def Filter(self, wavelength, data):
        out = []
        for enum, var in enumerate(wavelength):
            Sigma = self.Resolution(var)
            #print var, Sigma
            new = filt.gaussian_filter1d(data, sigma=Sigma)
            new = new[enum]
            out.append(new)
        return out

    def smoothCurve(self, **kwargs):
        photon = self.data
        self.wgth, photon = tb.clean(self.wgth, photon)
        photon = self.Filter(self.wgth, photon)

        if (self.plot == True):
            fig = pl.figure()
            pl.plot(self.wgth, photon, color='b',
                    label='Convolved with seeing(lambda)')
            pl.xlabel('wavelength (nm)')
            pl.legend()
            pl.show()
        return self.wgth, photon
