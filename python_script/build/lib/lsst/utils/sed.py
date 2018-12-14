#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
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
from . import toolbox as tb

# ===========================================================


class SED(Exception):
    pass
# ===========================================================


def MoffatFit(pixels, projected_footprint, A, mu, sigma):

    fit_g = fitting.LevMarLSQFitter()
    psf = fit_g(g_init, pixels, projected_footprint)
    start = psf.x_0 - 5 * psf.gamma
    end = psf.x_0 + 5 * psf.gamma
    integral = (integrate.quad(lambda pixels: psf(pixels), start, end))[0]
    print("psf_moffat = ", integral, psf.x_0, psf.gamma)

    ''' begin Control plot'''
    pl.plot(pixels, psf(pixels), label='Fitted data')
    pl.yscale('log')
    pl.ylim(1., 1E6)
    pl.plot(projected_footprint)
    pl.show()
    ''' End control plot'''
    return


'''Match target name with its calibrated spectrum'''


class SED(object):
    def __init__(self, object_name, convolve=2., ronchi=200., pix2wgth=2.1,
                 method='gauss',
                 plot=False,
                 seeing_at_wght=[2.3e-06, -3.9e-03, 3.6], **kwargs):
        self.pix2wgth_coef = pix2wgth
        self.plot = plot
        self.calspec = False
        if ronchi is None:
            self.ronchi = 200.
        else:
            self.ronchi = float(ronchi)
        self.object_name = object_name
        self.Init()
        self.method = method
        self.seeing_at_wght = seeing_at_wght
        self.convolve = convolve/self.step # for constant seeing approx.
        self.Calspec()

    def Init(self):
        ''' this is for .9-m telescope '''
        self.step = 0.25      # SED wght step
        f_D = 13.5      # telescope f number
        l = 5.3*10**7 # distance between grating and focal plan (filter wheel 2)
        self.a = 1.*10**6 / self.ronchi  # default : 200 lines per mm gratings
        '''I tweak the known value to get a better modeling of resolution(wght)'''
        #self.a             = 8.*10**5 / self.ronchi  # default : 200 lines per mm gratings
        logging.info("ronchi period :"+str(self.a)+' nm')
        self.l_f_D = l / f_D
        self.pixel = 24*10**3

    def CalibSpec(self):
        ''' Load environment variables'''
        tb.ProdInstall()
        stellar_path = os.environ['SED_PREFIX']
        if (self.object_name == 'HD205905'):
            file_name = 'hd205905_stis_003.fits'
        elif (self.object_name.lower().find('hd14943') >= 0):
            file_name = 'hd14943_stis_003.fits'
        elif (self.object_name == 'HD200654'):
            file_name = 'hd200654_stis_003.fits'
        elif (self.object_name == 'HD185975'):
            file_name = 'hd185975_stis_003.fits'
        elif ((self.object_name == 'MUCOL') or (self.object_name == 'MuCol')):
            file_name = 'mucol_stis_003.fits'
        elif (self.object_name == 'HR7950'):
            file_name = 'hr7950.dat'
        elif (self.object_name == 'HR9087'):
            file_name = 'hr9087.dat'
        elif ((self.object_name == 'HD108344_disp') or (self.object_name == 'HD108344')):
            file_name = 'pickles_uk_2.ascii'
        elif ((self.object_name.lower()).find('lamlep') >= 0):
            file_name = 'lamlep_stis_004.fits'
        elif ((self.object_name.lower()).find('ksi02cet') >= 0):
            file_name = 'ksi2ceti_stis_004.fits'
        else:
            file_name = None
            logging.debug(str("No calibrated spectrum found for "+self.object_name))
            return file_name

        logging.info("Loading Calspec : "
                     + os.path.join(stellar_path, file_name))
        return os.path.join(stellar_path, file_name)

    '''Extract infos from Calspec fits file'''

    def Calspec(self):
        self.calib_spec = self.CalibSpec()
        if (self.calib_spec == None):
            self.calspec = False
            return
        else:
            self.calspec = True
            "All data are in Flam unit. Convert Flam into Photons"
            if (self.calib_spec.split('.')[-1] == 'fits'):
                fits = pf.open(self.calib_spec)
                tbdata = fits[1].data                  # the table is in extension 1
                cols = tbdata.columns                # names of the columns
                logging.debug(repr(fits[0].header))      # show fits header
                logging.debug(cols)                      # show table header
                tbdata = tbdata[tbdata.field(0) < 11000] # Select lambda<1100nm
                self.wavelength = tbdata.field('WAVELENGTH')/10        # Convert in nm
                self.flux = tbdata.field('FLUX') * self.wavelength
                self.staterror = tbdata.field('STATERROR') \
                    * self.wavelength * self.wavelength
                self.syserror = tbdata.field('SYSERROR')\
                    * self.wavelength * self.wavelength
            elif ((self.calib_spec.split('.')[-1] == 'dat')
                  or (self.calib_spec.split('.')[-1] == 'ascii')):
                spectrum = np.loadtxt(self.calib_spec)
                self.wavelength = spectrum[:, 0]/10
                self.flux = spectrum[:, 1] * self.wavelength
            """Rescale on 0.25 nm steps from 300 to 1100"""
            steps = np.arange(300, 1100.25, self.step)
            self.flux = interp.griddata(self.wavelength, self.flux, steps)
            self.wavelength = steps
            return

    ''' Geometric distorsion as a function of angle (wavelength)'''
    ''' due to defocus of a non collimated beam                 '''

    def Resolution(self, item, method='local'):
        delta = self.l_f_D * (1/(np.sqrt(1-(item/self.a)**2)) - 1)
        delta = delta / self.pixel # in fraction of pixel
        delta = delta * self.pix2wgth_coef
        focus = np.polyval(self.seeing_at_wght, item)
        focus = focus * self.pix2wgth_coef / self.step # convert in nm, then to the SED bin size

        '''if   : seeing taken from poly2 fit of psf_sigma array in raw_spectrum.list'''
        '''else : seeing taken from rawspectrum.list keyword (a constant)            '''
        if(method == 'local'):
            seeing = focus + delta
        elif(method == 'global'):
            seeing = self.convolve + delta
        return seeing   # in nm

    def Filter(self, wavelength, data):
        out = []
        for enum, var in enumerate(wavelength):
            if(self.method == 'gauss'):
                Sigma = self.Resolution(var)
                #print var, Sigma
                new = filt.gaussian_filter1d(data, sigma=Sigma)
            if(self.method == 'voigt'):
                ### must implement here the determination of fwhm_g and fwhm_l
                ### from the parameters of their fit as a function of wght.
                ### Not done yet because both parameters look very noisy,
                ### I would like to try a weight map.
                voigt = models.Voigt1D(amplitude_L=1., x_0=var, fwhm_G=fwhm_g,
                                       fwhm_L=fwhm_l)
                new = filt.generic_filter1d(data, voigt)
            new = new[enum]
            out.append(new)
        return out

    def smoothedSed(self, eval_above_telescope=False, **kwargs):
        photon = self.flux
        self.wavelength, photon = tb.clean(self.wavelength, photon)
        photon = self.Filter(self.wavelength, photon)

        if (self.plot == True):
            fig = pl.figure()
            norm = max(photon)
            phot = photon / norm
            pl.plot(self.wavelength, phot, color='b',
                    label='SED convolved with seeing(lambda)')
            pl.xlabel('wavelength (nm)')
            pl.ylabel('renorm (just for the plot)')
            pl.legend()
            pl.show()
        return self.wavelength, photon

    def SedModel(self, telescope, eval_above_telescope=False, **kwargs):
        wg, qe = tb.telescope_T(telescope)
        w_atm, T_atm = tb.AtmoT()
        photon = self.flux
        if (eval_above_telescope is True):
            qe = np.ones(len(self.wavelength))
        else:
            qe = interp.griddata(wg, qe, self.wavelength)
        T_atm = interp.griddata(w_atm, T_atm, self.wavelength)
        photon = photon * T_atm * qe
        # must clean list from nan before 1d filter.
        self.wavelength, photon = tb.clean(self.wavelength, photon)
        # for constant resolution :
        # photon       = filt.gaussian_filter1d(photon, sigma= self.convolve)
        photon = self.Filter(self.wavelength, photon)
        #clean        = np.array([self.wavelength, photon]).transpose()
        #clean        = (clean[~np.isnan(clean).any(1)]).transpose()
        #photon       = clean[1,:]
        norm = max(photon)
        photon /= norm
        #return clean[0,:], photon
        fig = pl.figure()
        pl.plot(self.wavelength, photon, color='b')
        pl.xlabel('wavelength (nm)')
        pl.ylabel('SED* T_atm * Tel. convolved with seeing')
        pl.title('HD14943 - forward model')
        pl.show()
        pl.legend()
        fig.savefig('model.pdf')
        return self.wavelength, photon

    def SedAboveAtmo(self):
        photon = self.flux
        photon = self.Filter(self.wavelength, photon)
        #photon       = filt.gaussian_filter1d(photon, sigma = self.convolve)#constant
        clean = np.array([self.wavelength, photon]).transpose()
        clean = (clean[~np.isnan(clean).any(1)]).transpose()
        return clean[0, :], clean[1, :]

    def QEMirror(self):
        wg, qe = tb.ccd_qe()
        wg2, R = tb.MirrorR()
        R = interp.griddata(wg2, R, wg)
        qeR = qe * R * R
        return wg, qeR

    '''Estimate the atmospheric transmission from dividing  '''
    '''Obs by T_tel by SED of reference star'''

    def EstimateAtmoT(self, wgth, signal, mkcalib=True, **kwargs):
        sedwg, sedflux = self.SedAboveAtmo()
        if mkcalib:
            wg, qeMr = self.QEMirror()
            qeMr = interp.griddata(wg, qeMr, wgth)
            signal = signal / qeMr
        sedflux = interp.griddata(sedwg, sedflux, wgth)
        atmo = signal / sedflux
        clean = np.array([wgth, atmo]).transpose()
        clean = (clean[~np.isnan(clean).any(1)]).transpose()
        return clean[0, :], clean[1, :]
