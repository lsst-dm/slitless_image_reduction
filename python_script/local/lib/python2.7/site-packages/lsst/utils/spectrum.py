#!/usr/bin/env python

import telinst as instru
import os
import sys
import astropy.io.fits as pf
import numpy as np
from croaks import NTuple
import toolbox as tb
import reduceimage as ri
import extraction as ex
import flatfield as fd
import dispersion as disp
import isr
import logging
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def frame(fits, x_start, x_end, y_start, y_end):
    return frame


# ===========================================================
class spectrum(Exception):
    pass
# ===========================================================


''' A class to extract a spectrum from a fits image'''
''' Warning : x and y from usual ds9 are exchanged due to python convention !!!'''


class spectrum(object):
    def __init__(self, object_name, x_object, y_object,
                 image,
                 x_start, x_end, y_start, y_end,
                 mask=None,
                 dispersion=False,
                 calibrate=False,
                 orientation='x',
                 plot=True,
                 Map=False,
                 out_dir='./',
                 seeing=None,
                 sx=None,
                 sy=None,
                 isr=False,
                 cosmic=False,
                 order='m+1',
                 method=1,
                 offset=100, **kwargs):

        self.calibrate = calibrate
        self.orientation = orientation
        self.dispersion = dispersion
        self.out_dir = out_dir
        self.plot = plot
        self.Map = Map
        self.isr = isr
        self.cosmic = cosmic
        self.method = method
        self.seeing = seeing
        self.sx = sx   # median sigma x of clump of stars
        self.sy = sy   # median sigma y of clump of stars
        self.object_name = object_name
        self.order = order
        self.offset = offset
        self.x_start, self.x_end, self.y_start, self.y_end\
            = int(x_start), int(x_end), int(y_start), int(y_end)
        self.x_object, self.y_object = x_object, y_object
        self.image, self.mask = image, mask
        self.LoadVar()

    ''' Load environment and prod variables'''

    def LoadVar(self):
        tb.ProdInstall()
        self.monitoring_photodiode_file =\
            os.path.join(os.environ['INSTRU_PREFIX'], 'S2281.tuple')
        self.jd = tb.JD(self.image)
        self.run = tb.Run(self.jd)
        self.filters = tb.Filters(self.image)
        if self.dispersion is True:
            self.dispers = disp.dispersion(self.run, self.filters, offset=self.offset)
        self.aper_flux = None
        self.calib_profile = None
        self.pixel = None
        self.wavelength = None
        if (self.order == 'm+1'):
            self.flip = False
        elif(self.order == 'm-1'):
            self.flip = True

    '''The actual extraction of the spectrum'''

    def TargetSpectrum(self, Direction='y', **kwargs):
        logging.info('Opening of image : ' + str(self.image))
        inst = instru.telinst(self.image, verbose='')
        image_data = inst.Image(self.image)
        head = inst.header
        #name       = head['OBJECT'] + head['RECID']
        footprint = tb.Flip(image_data[self.y_start: self.y_end,
                                       self.x_start: self.x_end], flip=self.flip)  # Select the footprint

        '''Loading frame from segmentation.fits'''
        if(self.mask):
            logging.info('Extracting frame from : ' + str(self.mask))
            seg = instru.telinst(self.mask, verbose='')
            m = seg.Image(self.mask)
            footprint_mask = tb.Flip(m[self.y_start: self.y_end,
                                       self.x_start: self.x_end],
                                     flip=self.flip)  # Flip footprint mask if m-1 order

        '''Subtract to the footprint a periodic electronic patterns in parallel direction'''
        profile = []
        if((self.mask) and (self.isr == True)):
            logging.info("Doing ISR")
            Isr = isr.isr(self)
            profile = Isr.Mask(inst, Direction)
            if (len(profile) != 0): # subtract periodic background
                footprint = tb.SubtractProfile(footprint, profile)

        '''Determine and write map of defects   '''
        '''if True, Correct for defective pixels'''
        if self.dispersion is True:
            "initialize flatfielding : needed both for DefectMap and calibrate"
            flat = fd.flatfield(self, footprint, head, self.dispers)
            if ((self.Map is True) and (flat.masterflat_list != None)):
                footprint = flat.divideByMap(footprint)
            else:
                logging.info("No flatfielding by map of defects")

        ''' Extract a map of cosmic in footprint'''
        ''' Write footprint before removing cosmics in case of other than aperture phot'''
        ''' Cosmics in footprint are also replaced with median of surrounding pixels'''
        if((self.mask) and (self.cosmic == True)):
            mean, sigma = tb.BkgdStat(footprint, footprint_mask)
            cosmics = ri.filters(plot=self.plot)
            cosmicimage = cosmics.Cosmics2(footprint, self.seeing, mean, sigma)
            hdu = pf.PrimaryHDU(cosmicimage)
            name = str("cosmics"+self.order+".fits")
            hdu.writeto(os.path.join(self.out_dir, name), overwrite=True)
            logging.info("Cosmics that are found are replaced with median of surrounding pixels ")

        hdu = pf.PrimaryHDU(footprint)
        name = str("footprint"+self.order+".fits")
        hdu.writeto(os.path.join(self.out_dir, name), overwrite=True)

        ''' Extract raw profile'''
        ''' Testing two versions:'''
        '''1st is constant aperture'''
        '''second used gaussian fit to determine center and width'''
        logging.info("Spectrum extraction method : "+str(self.method))

        get = ex.extract(self, footprint, plot=self.plot)

        if (self.method == 0):
            get.flux(mode='both')
            self.pixel, self.aper_flux = get.pixel, get.aper_flux

        if (self.method == 1):
            self.pixel, self.aper_flux = get.ExtractSpectrum()
            #self.pixel, self.aper_flux = get.aperture()

        if (self.method == 2):
            '''First pass to determine center from gaussian fit'''
            pixel, aper_flux, center, flux, sigma =\
                get.ExtractSpectrum2(mode='psf',
                                     plot=self.plot)

            '''The trace of the centroid of the spectrum'''
            traceC, traceS = get.centroidSpectrum(pixel, plot=True)
            '''Now, ready to re-run the spectrum extraction task for aperture, and with better guess for the psf fitting'''
            self.pixel, self.aper_flux, center, psf_flux, sigma =\
                get.ExtractSpectrum2(mode='both',
                                     plot=self.plot,
                                     position=traceC,
                                     aperture=traceS)

            '''could measure it here instead of in method 1 only '''
            get.psf_voigt_flux = np.zeros(len(self.pixel))
            get.psf_voigt_fwhmL = np.zeros(len(self.pixel))
            get.psf_voigt_fwhmG = np.zeros(len(self.pixel))

        self.pixel = self.FramePos2RefPos(self.pixel)

        '''Determine dispersion relation'''
        if self.dispersion is True:
            self.wavelength = self.dispers.Pixel2Wavelength(self.pixel)
            '''calibration using synthetic flat'''
            if self.calibrate is True:
                flat.Calibrate(footprint, head)
        else:
            self.wavelength = np.zeros(len(self.pixel))

        '''Dump raw spectrum'''
        self.DumpSpectrum(head,
                          ('pixel', 'w', 'aper_flux',
                           'psf_gauss_flux', 'psf_gauss_sigma', 'psf_gauss_mu',
                           'psf_moffat_flux', 'psf_moffat_x0', 'psf_moffat_gamma',
                           'integralGM',
                           'integralG',
                           'amplitude_0',
                           'x_0_0',
                           'gamma_0',
                           'alpha_0',
                           'amplitude_1',
                           'mean_1',
                           'stddev_1'),
                          (self.pixel, self.wavelength, np.array(self.aper_flux),
                           np.array(get.psf_gauss_flux),
                           np.array(get.psf_gauss_sigma),
                           np.array(get.psf_gauss_mu),
                           np.array(get.psf_moffat_flux),
                           np.array(get.psf_moffat_x0),
                           np.array(get.psf_moffat_gamma),
                           np.array(get.integralGM),
                           np.array(get.integralG),
                           np.array(get.amplitude_0),
                           np.array(get.x_0_0),
                           np.array(get.gamma_0),
                           np.array(get.alpha_0),
                           np.array(get.amplitude_1),
                           np.array(get.mean_1),
                           np.array(get.stddev_1)),
                          str("rawspectrum"+self.order+".list"))

        return

    def DumpSpectrum(self, head, names, list, file):
        list = zip(*list)
        info = np.rec.fromrecords([i for i in list], names=names)

        info = info.view(NTuple)
        if (self.filters.find('RONCHI400') >= 0):
            info.keys['RONCHI'] = 400
        else:
            info.keys['RONCHI'] = 200
        info.keys['IMG_NAME'] = str(head.get('IMG_NAME')+' ' + self.out_dir)
        info.keys['AIRMASS'] = head.get('AIRMASS')
        info.keys['OUTHUM'] = head.get('OUTHUM')
        info.keys['OUTPRESS'] = head.get('OUTPRESS')
        info.keys['WNDSPEED'] = head.get('WNDSPEED')
        info.keys['EXPTIME'] = head.get('EXPTIME')
        info.keys['HA'] = head.get('HA')
        info.keys['LATITUDE'] = '-30 10 07.90'
        info.keys['DEC'] = head.get('DEC')
        info.keys['ZENITH'] = head.get('ZD')
        info.keys['FILTERS'] = head.get('FILTERS')
        info.keys['X'] = self.x_object
        info.keys['Y'] = self.y_object
        info.keys['JD'] = self.jd
        info.keys['DEFECT_MAP'] = self.Map
        info.keys['DISPERSION'] = self.dispersion
        if self.dispersion is True:
            info.keys['PIX2WGTH'] = self.dispers.pix2wgth_coef[0]
        else:
            info.keys['PIX2WGTH'] = 0
        info.keys['FLATFIELDED'] = self.calibrate
        info.keys['ORIENTATION'] = self.orientation
        info.keys['SEEING'] = self.seeing
        info.keys['SX'] = self.sx
        info.keys['SY'] = self.sy
        info.keys['OBJECT'] = (self.object_name).lower()
        info.keys['INFO1'] = "pixel # start from direct object"
        info.totxt(os.path.join(self.out_dir, file))
        return

    ''' Position relative to direct image -> Position relative to footprint frame'''

    def RefPos2FramePos(self, pixel):
        pixel = np.asarray(pixel, dtype=np.float64)
        return pixel - self.offset

    ''' The other way around'''

    def FramePos2RefPos(self, pixel):
        pixel = np.asarray(pixel, dtype=np.float64)
        return pixel + self.offset
