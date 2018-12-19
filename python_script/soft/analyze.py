#!/usr/bin/env python
'''
return wght calibrated spectrum + forward model + EW measurements

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

from builtins import zip
from builtins import str
from builtins import object
import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import scipy.interpolate as interp
import sed as Sed
import lines as Lines
import argparse
import logging


'''-------------------------------------'''


class Modeling(object):
    def __init__(self, target, ronchi,
                 convolve=2.,
                 directory='./',
                 plot=False,
                 pix2wgth=2.1,
                 telescope='open',
                 seeing_at_wght=[2.3e-06, -3.9e-03, 3.6],
                 eval_above_telescope=False, **kwargs):
        self.target = target
        self.ronchi = ronchi
        self.directory = directory
        self.plot = plot
        self.pix2wgth_coef = pix2wgth
        self.convolve = convolve
        self.eval_above_telescope = eval_above_telescope
        self.telescope = telescope
        self.Init(seeing_at_wght)

    def Init(self, seeing_at_wght):
        self.sed = Sed.SED(self.target, self.convolve,
                           pix2wgth=self.pix2wgth_coef,
                           ronchi=self.ronchi,
                           seeing_at_wght=seeing_at_wght)
        if self.sed.calspec:
            self.sed_wavelength, self.sed_flux = self.sed.SedModel(self.telescope)
            if (self.eval_above_telescope is True):
                self.sed_wavelength_noqe, self.sed_flux_noqe\
                    = self.sed.SedModel(eval_above_telescope=True)

    '''Inherited from SED class'''

    def Resolution(self, wght):
        return self.sed.Resolution(wght)

    def ForwardModeling(self, wavelength, signal, name='', **kwargs):
        fig = pl.figure()
        signal = np.array(signal) # to convert in unit of energy
        norm = max(signal[~np.isnan(signal)])
        signal /= norm
        pl.plot(wavelength, signal, color='b', label='Obs ctio 0.9-m')
        if self.sed.calspec:
            #sed_wavelength, sed_flux = self.sed.SedModel()
            pl.plot(self.sed_wavelength, self.sed_flux, '-k', label='SED forward modeling')

        pl.xlabel('wavelength (nm)')
        pl.ylabel('flux (norm@max)')
        pl.title(str(self.target+' '+name))
        ymin = min(min(self.sed_flux), np.nanmin(signal))
        ymax = max(max(self.sed_flux), np.nanmax(signal))
        xmin = min(np.nanmin(self.sed_wavelength), np.nanmin(wavelength))
        xmax = max(np.nanmax(self.sed_wavelength), np.nanmax(wavelength))
        xmax = min(xmax, 1100)
        pl.xlim(xmin, xmax)
        pl.ylim(ymin, ymax)
        pl.legend()
        fig.savefig(os.path.join(self.directory, self.target+name+'.pdf'))

    def Atmo(self, wavelength, signal, mkcalib=True, name='', **kwargs):
        if self.sed.calspec:
            fig = pl.figure()
            signal = np.array(signal)
            wg, atmo = self.sed.EstimateAtmoT(wavelength, signal, mkcalib=True)
            tb.DumpTuple(['wg', 'atmo'],
                         [wg, atmo],
                         os.path.join(self.directory, str('Tatmo'+name+'.list')))
            if self.plot is True:
                pl.plot(wg, atmo, '-k', label=str('Tatmo'+name))
                pl.xlabel('wavelength (nm)')
                pl.ylabel('Atmospheric transmission (norm@max)')
                ymin = np.nanmin(atmo)
                ymax = np.nanmax(atmo)
                xmin = np.nanmin(wg)
                xmax = np.nanmax(wg)
                xmax = min(xmax, 1100)
                pl.xlim(xmin, xmax)
                pl.ylim(ymin, ymax)
                pl.legend()
        return


'''-------------------------------------'''


def func(x, a, b):
    return a * np.arcsin(x/b)


'''fit using features position found in both the obs. and the template'''


def Match(model, obs):
    interval = 6.
    Lm = []
    Lo = []
    for line in model:
        for i in obs:
            if((i >= line-interval) and (i <= line+interval)):
                Lm.append(line)
                Lo.append(i)
    return Lm, Lo


'''Return the coefficients to adjust the wght solution on the template'''
'''Aslo return the list of features detected on both the template and the obs'''


def Fit(model_file, obs_file, plot):
    obs_values = tb.readtxt(obs_file, ['wg_lines', 'flux_lines', 'sflux_lines'])
    model_values = tb.readtxt(model_file, ['wg_lines', 'flux_lines', 'sflux_lines'])
    Lm, Lo = Match(model_values[0], obs_values[0])
    if plot is True:
        fig, ax1 = pl.subplots()
        pl.xlabel('obs')
        pl.ylabel('template')
        ax1.plot(Lo, Lm, 'r^')

    degree = 2
    if (len(Lm) <= degree):
        degree = int(len(Lm)-1)
        print("Changed degree to : ", degree)

    fit_coef = np.polyfit(Lo, Lm, deg=degree)
    points = np.linspace(100, 1500, 100)
    dispersion = np.polyval(fit_coef, Lo)
    dispersion2 = np.polyval(fit_coef, points)
    if plot is True:
        ax1.plot(points, dispersion2, color='k')
        fig, ax1 = pl.subplots()
        ax1.plot(Lo, Lm-dispersion, 'r^')
        pl.axhline(y=0., linewidth=1, color='k')
        pl.xlabel('template')
        pl.ylabel('residuals')
    print('rms residuals = ', np.std(Lm-dispersion))
    return fit_coef, Lm


'''Extraction of EW from automated feature location'''


def Extraction(name, wavelength, flux):
    data = np.array([wavelength, flux]).transpose()
    li = Lines.Lines(data)
    wg_lines, flux_lines, sflux_lines = li.ExtractLines()
    li.Plot(write=True,
            name=str(name+'.pdf'))
    tb.DumpTuple(['wg_lines', 'flux_lines', 'sflux_lines'],
                 [wg_lines, flux_lines, sflux_lines],
                 str(name+'.list'))


'''Extraction of EW from template input'''


def EW(wavelength, flux, positions, sigmas):
    ew = []
    segment = []
    data = np.array([wavelength, flux]).transpose()
    for i, j in zip(positions, sigmas):
        cut = data[(data[:, 0] >= (i-j)) & (data[:, 0] <= (i+j))]
        linear_coef = np.polyfit([cut[0, 0], cut[-1, 0]],
                                 [cut[0, 1], cut[-1, 1]], deg=1)
        bkgd = np.polyval(linear_coef, cut[:, 0])
        segment.append([cut[:, 0], bkgd])
        equ_w = np.trapz(1 - cut[:, 1] / bkgd)
        print('pos, width, EW : ', i, 2*j, equ_w)
        ew.append(equ_w)
    return ew, segment


def Plot(wavelength, flux, segment, name='lines2.pdf', **kwargs):
    fig = pl.figure()
    pl.plot(wavelength, flux, color='black', label='Obs')
    for i in segment:
        pl.plot(i[0], i[1], color='blue', linewidth=2)
    pl.title('Lines extraction')
    pl.xlabel('Wavelength (nm)')
    pl.ylabel('Spectrum Profile (arbitrary units)')
    pl.legend(loc='upper right')
    fig.savefig(name)
    return


'''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
'''to dispersion direction (return poly2 fit parameters)                         '''


def SeeingAtWght(w, psf_sigma):
    clean = np.array([w, psf_sigma]).transpose()
    clean = clean[(clean[:, 0] > 400) & (clean[:, 0] < 900)]
    seeing_at_wght = np.polyfit(clean[:, 0], clean[:, 1], deg=2)
    return seeing_at_wght


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "process slitless images"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="analyze spectrum")
    parser.add_argument('-p', "--plot",
                        help="show control plots",
                        action='store_true')
    parser.add_argument('-e', "--eval",
                        help="evaluate SED above telescope",
                        action='store_true')
    parser.add_argument('-v', "--verbose",
                        help="verbose",
                        action='store_true')
    parser.add_argument('-f', "--file", type=str,
                        help="name of spectrum file",
                        default='rawspectrum.list')
    parser.add_argument('-n', "--name", type=str,
                        help="output files identifiers",
                        default='+1')
    parser.add_argument('-i', "--imgs", nargs='+', type=str,
                        help="list of spectrum to be analyzed",
                        default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    # Run the atmo after the recalibration of the disperison relation.

    # Dump a spectrum.list with the old and new dispersion, and a model column.
    # Also dump a pdf comparing the template with the recalibrated observation.
    # From the recalib. dump new EW estimate from fix edges positions.

    args = grabargs()
    plot = args.plot
    verbose = args.verbose
    input_reps = args.imgs
    file = args.file
    name = args.name
    eval_above_telescope = args.eval
    if(verbose):
        Level = logging.getLogger().setLevel(logging.DEBUG)
        logging.debug('DEBUG mode')
    else:
        Level = logging.getLogger().setLevel(logging.INFO)
        logging.info('INFO mode')

    tb.ProdInstall()

    template_path = os.environ['TEMPLATE_PREFIX']
    if not os.path.exists(template_path):
        os.makedirs(template_path)

    for input_rep in input_reps:
        logging.info('\n' + 'Analyzing : ' + input_rep)

        '''observed spectrum'''
        param_file = os.path.join(input_rep, file)
        if not os.path.isfile(param_file):
            continue
        array, keys = tb.readlist(param_file,
                                  ['pixel', 'w',
                                   'aper_flux',
                                   'psf_gauss_flux', 'psf_gauss_sigma', 'psf_gauss_mu'])
        target = tb.STDname(keys.get('OBJECT'))
        seeing = keys.get('SEEING') # This is currently not a good measurement
        pix2wgth = keys.get('PIX2WGTH')
        ronchi = keys.get('RONCHI')
        array = np.array(array).transpose()
        array = array[(array[:, 1] <= 1100)] # Need a better model for the fitting
        pixel = array[:, 0]
        w = array[:, 1]
        aper_flux = array[:, 2]
        psf_gauss_flux = array[:, 3]
        psf_gauss_sigma = array[:, 4] # now use for a better seeing determination as a function of wght
        psf_gauss_mu = array[:, 5]

        ''' determine telescope filter set up'''
        filters = keys.get('FILTERS')
        print(filters)
        if (str(filters).find('RG715') >= 0):
            print('RG715 filter in place')
            tel_t = 'rg715'
        else:
            print('Telescope set up open')
            tel_t = 'open'

        '''Select the best flux estimator'''
        #best_flux = psf_gauss_flux #psf_voigt_flux
        best_flux = aper_flux #psf_voigt_flux

        '''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
        '''to dispersion direction (return poly2 fit parameters)                         '''
        seeing_at_wght = SeeingAtWght(w, psf_gauss_sigma)

        '''Extraction of absorption lines in observation'''
        obs_file = os.path.join(input_rep, str('obs_lines'+name))
        Extraction(obs_file, w, best_flux)

        ''' Load CalSpec reference spectrum based on target name   '''
        ''' Accounting for atmospheric transmission and resolution '''
        ''' Atmospheric transmission from forward modeling         '''

        convolve = pix2wgth * seeing # in nm # constant approximation, deprecated

        mod = Modeling(target, ronchi,
                       convolve=convolve,
                       directory=input_rep,
                       plot=True,
                       pix2wgth=pix2wgth,
                       seeing_at_wght=seeing_at_wght,
                       eval_above_telescope=eval_above_telescope)

        mod.ForwardModeling(w, best_flux, name='_1')
        mod.Atmo(w, best_flux, mkcalib=True, name='_1')

        '''Extraction of absoption lines in template '''
        model_file = os.path.join(template_path, str(target+'_lines'))
        if not os.path.isfile(str(model_file+'.list')):
            if mod.sed_wavelength is not None:
                Extraction(model_file, mod.sed_wavelength, mod.sed_flux)

        ''' Refit the dispersion relation by matching absorption features '''
        ''' between the observation and a template                        '''
        fit_coef, ew_position = Fit(str(model_file+'.list'), str(obs_file+'.list'), plot)

        ''' Shifting the initial wght calibration '''
        wght_recal = np.polyval(fit_coef, w)
        if plot is True:
            fig, ax1 = pl.subplots()
            ax1.plot(w, best_flux, 'r', label='before')
            ax1.plot(wght_recal, best_flux, 'k', label='after')
            pl.xlabel('wght')
            pl.ylabel('best_flux')

        '''SED above telescope, to estimate telescope throughput'''
        sed_flux_noqe = np.zeros(len(wght_recal))
        if (eval_above_telescope is True):
            sed_flux_noqe =\
                interp.griddata(mod.sed_wavelength_noqe, mod.sed_flux_noqe, wght_recal)

        '''Write analysis results in file'''
        sed_flux = interp.griddata(mod.sed_wavelength, mod.sed_flux, wght_recal)
        outfile = str('spectrum'+name+'.list')
        tb.DumpFile(keys, ('pixel', 'wght', 'wght_recal', 'flux', 'sed_flux', 'sed_flux_noqe'),
                    (pixel, w, wght_recal, best_flux, sed_flux, sed_flux_noqe),
                    os.path.join(input_rep, outfile))

        '''Measurement of the EW based on the position of the lines'''
        '''From the template +/- the resolution at this position '''
        sigmas = []
        segments = []
        for i in ew_position:
            sigmas.append(mod.Resolution(i))
        ew_flux, segments = EW(wght_recal, best_flux, ew_position, sigmas)
        outfile = str('ew_extraction'+name+'.pdf')
        Plot(wght_recal, best_flux, segments,
             name=os.path.join(input_rep, outfile))
        outfile = str('ew_extraction'+name+'.list')
        tb.DumpTuple(['wg_lines', 'flux_lines', 'sflux_lines'],
                     [ew_position, ew_flux, np.zeros(len(ew_position))],
                     os.path.join(input_rep, outfile))

        '''Replot after having rescaled the wght solution'''
        mod.ForwardModeling(wght_recal, best_flux, name='_2')
        mod.Atmo(wght_recal, best_flux, mkcalib=True, name='_2')

        if plot:
            pl.show()
        pl.clf()
