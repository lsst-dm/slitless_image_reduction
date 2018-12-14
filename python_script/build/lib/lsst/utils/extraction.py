#!/usr/bin/env python

from __future__ import print_function
from builtins import zip
from builtins import next
from builtins import str
from builtins import range
from builtins import object
import os
import sys
import matplotlib.cm as cm
import numpy as np
from scipy import integrate
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import astropy.io.fits as pf
import logging

#http://docs.astropy.org/en/stable/modeling/


# Define model function to be used to fit to the data above:
def Gauss1D(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


'''subtract mean of upper and lower bkgd of slice'''


def subtractBkgd(slice, width, bkgd_size):
    sub_footprint = slice[bkgd_size:(-1-bkgd_size+1)]
    bkgd1 = slice[:bkgd_size]        # ok
    bkgd2 = slice[width-bkgd_size:] # ok
    bkgd1 = np.median(bkgd1)
    bkgd2 = np.median(bkgd2)
    bkgd = np.mean([bkgd1, bkgd2])
    sub_footprint = sub_footprint - bkgd
    return sub_footprint


def MoffatFit(pixels, projected_footprint, A, mu, sigma, residuals, i, plot):
    g_init = models.Moffat1D(amplitude=A, x_0=mu, gamma=sigma)
    fit_g = fitting.LevMarLSQFitter()
    psf = fit_g(g_init, pixels, projected_footprint)

    start = psf.x_0 - 5 * psf.gamma
    end = psf.x_0 + 5 * psf.gamma
    integral = (integrate.quad(lambda pixels: psf(pixels), start, end))[0]
    ''' begin Control plot'''
    if((not i % 10)and(i < 400) and (plot == True)):
        pl.plot(pixels, psf(pixels), label='Moffat')
        pl.yscale('log')
        pl.ylim(1., 1E6)
        pl.plot(pixels, projected_footprint)
        pl.legend()
        pl.show()
        ''' End control plot'''
        '''Filling residuals'''
    #residuals[:, i] = psf(pixels)-projected_footprint
    return integral, psf.x_0.value, psf.gamma.value, psf.alpha.value


def GaussMoffatFit(pixels, projected_footprint,
                   amplitude, x_0, gamma, alpha,
                   A, mu, sigma,
                   residuals, i, plot):
    g_init = models.Moffat1D(amplitude=amplitude, x_0=x_0, gamma=gamma, alpha=alpha)\
        + models.Gaussian1D(amplitude=A, mean=mu, stddev=sigma)
    g_init.stddev_1.min = 0.5
    g_init.stddev_1.max = 3.
    g_init.amplitude_0.min = 1.
    g_init.amplitude_0.max = A/10.
    g_init.gamma_0.min = 1.
    g_init.gamma_0.max = 2.
    g_init.alpha_0.min = 1.
    g_init.alpha_0.max = 2.

    fit_g = fitting.LevMarLSQFitter()
    psf = fit_g(g_init, pixels, projected_footprint)

    start = psf.x_0_0 - 10 * psf.stddev_1
    end = psf.x_0_0 + 10 * psf.stddev_1
    integralGM = (integrate.quad(lambda pixels: psf(pixels), start, end))[0]
    integralG = np.sqrt(2 * np.pi) * psf.stddev_1 * psf.amplitude_1
    ''' begin Control plot'''
    if((not i % 10)and(i < 400) and (plot == True)):
        pl.plot(pixels, psf(pixels), label='Gauss+Moffat')
        pl.yscale('log')
        pl.ylim(1., 1E6)
        pl.plot(pixels, projected_footprint)
        pl.legend()
        pl.show()
    ''' End control plot'''
    '''Filling residuals'''
    residuals[:, i] = psf(pixels)-projected_footprint
    return integralGM, integralG, psf.amplitude_0.value, psf.x_0_0.value, psf.gamma_0.value, psf.alpha_0.value, psf.amplitude_1.value, psf.mean_1.value, psf.stddev_1.value


def VoigtFit(pixels, projected_footprint, A, mu, sigma, residuals, i, plot):
    g_init = models.Voigt1D(amplitude_L=A, x_0=mu, fwhm_G=sigma*2., fwhm_L=sigma*3)
    #g_init.x_0.min = 0.
    #g_init.x_0.max = 80.
    #g_init.fwhm_L.min = 0.1
    #g_init.fwhm_G.min = 0.1
    #g_init.fwhm_L.max = 20.
    #g_init.fwhm_G.max = 20.
    fit_g = fitting.LevMarLSQFitter()
    psf = fit_g(g_init, pixels, projected_footprint)
    start = psf.x_0 - 5 * (psf.fwhm_L+psf.fwhm_G)/2
    end = psf.x_0 + 5 * (psf.fwhm_L+psf.fwhm_G)/2
    integral = (integrate.quad(lambda pixels: psf(pixels), start, end))[0]
    ''' begin Control plot'''
    if((not i % 10)and(i < 400) and (plot == True)):
        pl.plot(pixels, psf(pixels), label='Voigt')
        pl.yscale('log')
        pl.ylim(1., 1E6)
        pl.plot(pixels, projected_footprint)
        pl.legend()
        pl.show()
    ''' End control plot'''
    return


# ===========================================================
class extract(Exception):
    pass
# ===========================================================


class extract(object):
    def __init__(self, spectrum, footprint, plot=False, **kwargs):
        self.footprint = footprint
        self.plot = plot
        self.spectrum = spectrum
        if self.spectrum.orientation == 'y':
            self.footprint = self.footprint.T
        logging.info("Extraction orientation : "+str(spectrum.orientation))
        self.psf_gauss_flux = None
        self.psf_gauss_sigma = None
        self.psf_gauss_mu = None
        self.init()

    def init(self):
        self.length = len(self.footprint[0])
        self.psf_gauss_flux = np.zeros(self.length)
        self.psf_gauss_mu = np.zeros(self.length)
        self.psf_gauss_sigma = np.zeros(self.length)
        self.pixel = np.zeros(self.length)
        self.aper_flux = np.zeros(self.length)
        self.psf_moffat_flux = np.zeros(self.length) #this is integrated on the pixel range
        self.psf_moffat_x0 = np.zeros(self.length)
        self.psf_moffat_gamma = np.zeros(self.length)

        self.integralGM = np.zeros(self.length)
        self.integralG = np.zeros(self.length)
        self.amplitude_0 = np.zeros(self.length)
        self.x_0_0 = np.zeros(self.length)
        self.gamma_0 = np.zeros(self.length)
        self.alpha_0 = np.zeros(self.length)
        self.amplitude_1 = np.zeros(self.length)
        self.mean_1 = np.zeros(self.length)
        self.stddev_1 = np.zeros(self.length)
        return

    '''Just a constant width aperture - minus mean bkgd from edges'''

    def flux(self, mode='aperture', **kwargs):
        bkgd_size = 10
        pixel = []
        aper_flux = []
        width = len(self.footprint)

        ysize = width - 2*bkgd_size
        xsize = self.length
        residuals = np.zeros([ysize, xsize])

        if((mode == 'psf') or (mode == 'both')):
            sigma = 3.
            mu = width/2.
            amplitude = 1000.

        for i in range(self.length):
            self.pixel[i] = i
            projected_footprint = subtractBkgd(self.footprint[:, i], width, bkgd_size)

            if((mode == 'aperture') or (mode == 'both')):
                Sum = np.sum(projected_footprint)
                self.aper_flux[i] = Sum
            if(((mode == 'psf') or (mode == 'both')) and (Sum > 1000.)):
                pixels = np.arange(0+bkgd_size, (width-bkgd_size))
                if ((mu > .7*width) or (mu < 0.3*width)):
                    mu = width/2.
                if ((sigma > 20.) or (sigma < 0.1)):
                    sigma = 3.
                p0 = [amplitude, mu, sigma]   # (A, mu, sigma)

                try:
                    coeff, var_matrix = curve_fit(Gauss1D, pixels, projected_footprint,
                                                  p0=p0,
                                                  bounds=(0., [1E6, 100., 10.]))
                    self.psf_gauss_sigma[i] = coeff[2]
                    self.psf_gauss_mu[i] = coeff[1]
                    amplitude = coeff[0]
                    perr = np.sqrt(np.diag(var_matrix))
                    self.psf_gauss_flux[i] = np.sqrt(
                        2 * np.pi) * self.psf_gauss_sigma[i] * amplitude # Gaussian integral
                    mu = coeff[1]   #starting value for next fit
                    sigma = coeff[2]   #starting value for next fit

                    '''Fitting a Moffat profile'''
                    #if((not i % 10) and (i<400) and (amplitude >=100)):
                    self.psf_moffat_flux[i],\
                        self.psf_moffat_x0[i],\
                        self.psf_moffat_gamma[i],\
                        alpha\
                        = MoffatFit(pixels, projected_footprint, amplitude, mu, sigma, residuals, i, self.plot)

                    '''Fitting Moffat+Gauss'''
                    self.integralGM[i], self.integralG[i],\
                        self.amplitude_0[i], self.x_0_0[i], self.gamma_0[i],\
                        self.alpha_0[i], self.amplitude_1[i], self.mean_1[i], self.stddev_1[i]\
                        = GaussMoffatFit(pixels, projected_footprint,
                                         self.psf_moffat_flux[i],
                                         self.psf_moffat_x0[i],
                                         self.psf_moffat_gamma[i],
                                         alpha,
                                         amplitude, mu, sigma,
                                         residuals, i, self.plot)

                except RuntimeError:
                    print('RuntimeError - moving on')

                '''Filling in the residuals'''
                fit = Gauss1D(pixels, *coeff)
                #residuals[:, i] = fit-projected_footprint #currently at Moffat

                ''' begin Control plot'''
                if(self.plot == True):
                    if((not i % 10)and(i < 400)):
                        print('aper : ', i, self.aper_flux[i])
                        print(i, self.psf_gauss_flux[i], self.psf_gauss_sigma[i], self.psf_gauss_mu[i], amplitude)
                        #fit = Gauss1D(pixels, *coeff)
                        pl.xlabel('Spectrum spatial profile (pixel)')
                        pl.ylabel('Amplitude (ADU)')
                        pl.title('CTIO .9m - %s'%(self.spectrum.object_name))
                        pl.plot(pixels, fit, label='Gauss')
                        pl.yscale('log')
                        pl.ylim(1., 1E6)
                        pl.plot(pixels, projected_footprint)
                        pl.legend()
                        pl.show()
                        ''' End control plot'''

        if((mode == 'psf') or (mode == 'both')):
            hdu = pf.PrimaryHDU(residuals)
            name = "residuals.fits"
            hdu.writeto(os.path.join(self.spectrum.out_dir, name), overwrite=True)
            logging.info("writing map of residuals")

        return

    '''Just a constant width aperture - minus mean bkgd from edges'''

    def apertureOld(self):
        bkgd_size = 10
        pixel = []
        aper_flux = []
        width = len(self.footprint)
        end = len(self.footprint[0])

        for i in range(end):
            projected_footprint = self.footprint[bkgd_size:(-1-bkgd_size+1), i] #ok
            bkgd1 = self.footprint[:bkgd_size, i]        # ok
            bkgd2 = self.footprint[width-bkgd_size:, i] # ok
            bkgd1 = np.median(bkgd1)
            bkgd2 = np.median(bkgd2)
            bkgd = np.mean([bkgd1, bkgd2])
            projected_footprint = projected_footprint - bkgd
            Sum = np.sum(projected_footprint)
            pixel.append(i)
            aper_flux.append(Sum)

        self.fillZeros(len(pixel))
        return pixel, aper_flux

    '''  Extracting spectrum from array     '''

    def ExtractSpectrum(self, mode='aperture', **kwargs):
        pixel = []
        aper_flux = []
        psf_flux = []
        psf_sigma = []
        psf_mu = []

        #print 'study of a smooth longitudinal approximation to set the weight map'
        #for i in range(len(footprint[0])):
        #    projected_footprint = footprint[i , :]
        #    pl.plot(projected_footprint)
        #    pl.show()

        sigma = 2. # starting sgma**2 value for fit
        amplitude = 1.
        extend = len(self.footprint)
        mu = extend / 2.
        pixels = np.arange(0, extend)
        end = len(self.footprint[0])
        for i in range(end):
            projected_footprint = self.footprint[:, i]
            pixel.append(i)
            '''Aperture flux'''
            Sum = np.sum(projected_footprint)
            aper_flux.append(Sum)
            if(np.isnan(projected_footprint).any()):#because there are nan after flatfielding
                continue

            '''Gaussian PSF flux, sigma and center'''
            if ((mu > 80.) or (mu < 0.)):
                mu = 40.
            if ((sigma > 20.) or (sigma < 0.1)):
                sigma = 2.
            p0 = [amplitude, mu, sigma]   # (A, mu, sigma)
            while True:
                try:
                    coeff, var_matrix = curve_fit(Gauss1D, pixels, projected_footprint,
                                                  p0=p0,
                                                  bounds=(0., [1E6, 100., 10.]))
                    break
                except RuntimeError:
                    print("Minimization failed")
                    break

            sigma = coeff[2]
            mu = coeff[1]
            amplitude = coeff[0]
            perr = np.sqrt(np.diag(var_matrix))
            integral = np.sqrt(2 * np.pi) * sigma * amplitude # Gaussian integral
            psf_flux.append(integral)
            psf_sigma.append(sigma)
            psf_mu.append(mu)

            #logging.debug('pixel, aper, fluxpsf, s_psf, center = %s %s %s %s %s'\
            #              %(i, Sum, integral, sigma, mu))

            #'''Fitting a Moffat profile'''
            #MoffatFit(pixels, projected_footprint, amplitude, mu, sigma)

            '''Fitting a Voigt profile       '''
            '''With weights from gaussian fit'''
            weights = Gauss1D(pixels, *coeff)**2
            a, b, c = VoigtFit(pixels, projected_footprint, amplitude, mu, sigma)#, weights=weights)
            self.psf_voigt_fwhmL.append(b)
            self.psf_voigt_fwhmG.append(c)
            #logging.debug('fluxvoigt, fwhmL, fwhmG = %s %s %s'\
            #              %(a,b,c))

            ''' begin Control plot'''
            #if(logging.getLogger().getEffectiveLevel()==10):
            #    fit = Gauss1D(pixels, *coeff)
            #    pl.xlabel('Spectrum spatial profile (pixel)')
            #    pl.ylabel('Amplitude (ADU)')
            #    pl.title('CTIO .9m - %s'%(self.object_name))
            #    pl.plot(pixels, fit, label='Gauss - pix:'+str(i))
            #    pl.yscale('log')
            #    pl.ylim(1., 1E6)
            #    pl.plot(projected_footprint)
            #    pl.legend()
            #    pl.show()
            #    ''' End control plot'''

        self.psf_gauss_flux = psf_flux
        self.psf_gauss_sigma = psf_sigma
        self.psf_gauss_mu = psf_mu

        return pixel, aper_flux

    '''  Extracting spectrum from array    mode are either [aperture, psf, both] '''

    def ExtractSpectrum2(self,
                         mode='aperture',
                         plot=True,
                         position=None,
                         aperture=None, # Default half spatial span of spectrum (pixel)
                         **kwargs):
        pixel = []
        aper_flux = []
        psf_flux = []
        psf_sigma = []
        psf_mu = []
        PIX = []
        PROJ = []

        amplitude = 1.
        extend = len(self.footprint)
        mu = extend / 2.
        pixels = np.arange(0, extend)
        sigma = 2.# starting sgma**2 value for fit

        if aperture is not None:
            logging.info('Variable aperture ')
            aperture = np.array(aperture)*6.#aper=6sigma
        else:
            aperture = 10.
            logging.info('Fixed aperture : ' + str(aperture))

        if position is not None:
            position = np.array(position)
            min = (np.floor(position-aperture)).astype(int)
            max = (np.ceil(position+aperture)).astype(int)

            #for i, j,k in zip(position, min,max):
            #    print 'array ', i, j,k

        end = len(self.footprint[0])
        for i in range(end):
            if position is not None:
                if(min[i] < 0):
                    min[i] = 0
                projected_footprint = self.footprint[min[i]:max[i], i]
            else:
                projected_footprint = self.footprint[:, i]
            pixel.append(i)
            '''Aperture flux'''
            if((mode == 'aperture') or (mode == 'both')):
                Sum = np.sum(projected_footprint)
                aper_flux.append(Sum)

            if((mode == 'psf') or (mode == 'both')):
                '''Gaussian PSF flux, sigma and center'''
                #if ((mu > 80.) or (mu<=0.)):
                if position is None:
                    mu = 40.
                else:
                    mu = (max[i]-min[i]) / 2. + min[i]
                if ((sigma > 20.) or (sigma < 0.1)):
                    sigma = 2.
                p0 = [amplitude, mu, sigma]   # (A, mu, sigma)
                if position is not None:
                    if min[i] < 0:
                        min[i] = 0
                    if max[i] > extend:
                        max[i] = extend
                    pixels = np.arange(min[i], max[i])
                    #if(plot==True):
                    #    print 'pixel :', i, ' min/max ', min[i], '/'\
                    #    ,max[i], len(pixels), ' ?= ', len(projected_footprint)
                if (len(pixels) == len(projected_footprint)):
                    sigma = 0.
                    mu = 0.
                    integral = 0.
                    #if (i>70) and (i<300):
                    if (i > 120) and (i < 180):
                        PIX.append(pixels)
                        PROJ.append(projected_footprint)
                    cpixels, cprojected_footprint = cleanZeros(pixels, projected_footprint)
                    try:
                        coeff, var_matrix = curve_fit(Gauss1D, cpixels, cprojected_footprint,
                                                      p0=p0,
                                                      sigma=1.,#np.sqrt(np.abs(cprojected_footprint)),
                                                      bounds=([0., 10., .5], [1E6, 100., 10.]),
                                                      check_finite=False)
                        sigma = coeff[2]
                        mu = coeff[1]
                        amplitude = coeff[0]
                        perr = np.sqrt(np.diag(var_matrix))
                        integral = np.sqrt(2 * np.pi) * sigma * amplitude # Gaussian integral

                    except RuntimeError:
                        print('RuntimeError - moving on')
                        #break

                psf_flux.append(integral)
                psf_sigma.append(sigma)
                psf_mu.append(mu)
                ''' begin Control plot'''
                check = False
                if((check == True) and (mode == 'both') and (i > 70) and (i < 300)):
                    fit = Gauss1D(pixels, *coeff)
                    pl.xlabel('Spectrum spatial profile (pixel)')
                    pl.ylabel('Amplitude (ADU)')
                    pl.title('CTIO .9m - %s'%(self.object_name))
                    pl.plot(pixels, fit, label='Gauss ')
                    pl.yscale('log')
                    pl.ylim(1., 1E6)
                    pl.plot(pixels, projected_footprint)
                    pl.legend()
                    pl.show()
                    ''' End control plot'''

        if (mode == 'psf'):
            showProfil(PIX, PROJ)

        print('mode ', mode, ' length ', len(pixel))
        len(aper_flux), len(psf_mu), len(psf_flux), len(psf_sigma)
        self.psf_gauss_flux = psf_flux
        self.psf_gauss_sigma = psf_sigma
        self.psf_gauss_mu = psf_mu

        return pixel, aper_flux, psf_mu, psf_flux, psf_sigma

    def centroidSpectrum(self, pixel, plot=False, **kwargs):
        clip = 1.5 # reject 2-sigma outliers
        center_clip = sigma_clip(self.psf_gauss_mu, sigma=clip)
        sub_pix = []
        sub_center = []
        sub_flux = []
        sub_sigma = []
        for i, j, k, l in zip(pixel, center_clip,
                              self.psf_gauss_flux, self.psf_gauss_sigma):
            if j != 0:
                sub_pix.append(i)
                sub_center.append(j)
                sub_flux.append(k)
                sub_sigma.append(l)

        f1 = np.polyfit(sub_pix, sub_sigma, 2, w=np.sqrt(sub_flux))
        f2 = np.polyfit(sub_pix, sub_center, 1, w=np.sqrt(sub_flux))
        seeing = np.poly1d(f1)
        p = np.poly1d(f2)
        if (plot == True):
            pl.clf()
            fig = pl.figure(1)
            pl.plot(pixel, self.psf_gauss_mu, 'bo', label='data')
            pl.plot(pixel, center_clip, 'r^', label='clip')
            pl.plot(pixel, p(pixel), color='blue', label='interp.')

            fig = pl.figure(2)
            pl.plot(pixel, self.psf_gauss_sigma, 'bo', label='data')
            pl.plot(sub_pix, sub_sigma, 'r^', label='clip')
            pl.plot(pixel, seeing(pixel), color='blue', label='interp.')
            pl.legend()
            pl.show()

        '''The trace of the centroid and sigma of the spectrum'''
        trace1 = p(pixel)
        trace2 = seeing(pixel)
        return trace1, trace2

    def fillZeros(self, length):
        zeros = [np.zeros(length)]
        self.psf_gauss_flux = np.zeros(length)
        self.psf_gauss_sigma = np.zeros(length)
        self.psf_gauss_mu = np.zeros(length)
        self.psf_voigt_flux = np.zeros(length)
        self.psf_voigt_fwhmL = np.zeros(length)
        self.psf_voigt_fwhmG = np.zeros(length)
        return

#    *************************    #


def cleanZeros(a, b):
    aa = []
    bb = []
    for i, j in zip(a, b):
        if j != 0.:
            aa.append(i)
            bb.append(j)
    return aa, bb


def showProfil(a, b):
    pl.figure()
    colors = iter(cm.rainbow(np.linspace(0, 1, len(a))))

    #pl.ylim(1., 1.2)
    for i, j in zip(a, b):
        norm = max(j)
        pl.plot(i, j/norm, color=next(colors), marker='x', linestyle='None')
        p0 = [1., 40., 2.]
        coeff, var_matrix = curve_fit(Gauss1D, i, j/norm, p0=p0)
        fit = Gauss1D(i, *coeff)
        #print i , fit
        pl.plot(i, fit, label='Gauss')
        pl.yscale('log')
    pl.xlabel('Spectrum spatial profile (pixel)')
    pl.ylabel('Amplitude (norm)')
    pl.show()
