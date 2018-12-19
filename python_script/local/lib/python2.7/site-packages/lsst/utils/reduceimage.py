#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from builtins import object
import os
import sys
import re
from . import telinst as instru
import numpy as np
import pylab as pl
from . import spectrum
import astropy.io.fits as pf
import logging
from . import toolbox as tb
import logging


'''----------------------------------------------------'''
''' A class to prepare raw fits for spectrum extraction'''
'''----------------------------------------------------'''


class Prepare(object):
    def __init__(self, Image,
                 object_name=None,
                 suffix='',
                 proddir=None,
                 **kwargs):
        self.img = Image
        self.object_name = object_name
        self.suffix = suffix
        self.proddir = proddir
        self.LoadVar()

    def LoadVar(self):
        ''' Load environment variables'''
        tb.ProdInstall(prodpath=self.proddir)
        prodpath = os.environ['PROD_PREFIX']
        self.instrument = instru.telinst(self.img)
        if self.object_name is None:
            self.object_name = str(self.instrument.header.get('OBJECT'))
            self.object_name = self.object_name.replace(' ', '_')
        self.hdr = (self.instrument.header).copy()
        out_rep_name = str(self.hdr['DATE']).replace('-', '_').replace(':', '_')+str(self.suffix)
        self.outdir = os.path.join(prodpath, self.object_name, out_rep_name)
        jd = tb.JD(self.img)
        self.run = tb.Run(jd)
        logging.info('JD :' + str(jd))
        logging.info('run : ' + str(self.run))
        return

    def Check(self):
        exptime = self.hdr['EXPTIME']
        if float(exptime) == 0.:
            logging.info('Image is not what you think. No processing done.')
            return True
        else:
            return False

    def TrimOverscanSubtract(self):
        '''overscan subtract, trim and write the image in PROD_PREFIX'''
        data = self.instrument.Image(self.img)
        outimg = self.instrument.OverscanSubtract_andTrim(data)
        self.hdr.add_comment("Image is trimmed")
        (filepath, self.filename) = os.path.split(self.img)
        self.hdr['IMG_NAME'] = (self.filename, 'initial image name')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.outname = os.path.join(self.outdir, self.filename)
        pf.writeto(self.outname, outimg, self.hdr, overwrite=True)
        return

    ''' Subtract master-dark -- scaled to same exposure time'''
    ''' could indicate the period and find dark or do nothing'''

    def Dark(self, dark):
        calibrated = os.path.join(self.outdir, 'calibrated.fits')
        coef = 1.
        found_dark = False
        if dark is True:
            if((self.run == 1) or (self.run == 2)):
                master_dark = os.path.join(os.environ['MASTER_PREFIX'],
                                           'master_dark_n1.fits')
                coef = self.hdr['EXPTIME'] / (pf.open(master_dark))[0].header['EXPTIME']
                found_dark = True
            if(self.run == 3):
                master_dark = os.path.join(os.environ['MASTER_PREFIX'],
                                           'masterbiasNov2016.fits')
                found_dark = True
            if(self.run == 5):
                master_dark = os.path.join(os.environ['MASTER_PREFIX'],
                                           'MasterBias_run5.fits')
                found_dark = True
            if(self.run == 7):
                master_dark = os.path.join(os.environ['MASTER_PREFIX'],
                                           'masterbiasOct2017.fits')
                found_dark = True
            if(found_dark == True):
                cmd = 'compute.py %s + -%s  %s %s'\
                      % (self.outname, coef, master_dark, calibrated)
                logging.info(cmd)
                os.system(cmd)
                os.remove(os.path.join(self.outdir, self.filename))
            else:
                print('no masterdark/bias found for this run')
                os.rename(self.outname, calibrated)
        else:
            os.rename(self.outname, calibrated)
        return

    '''divide by masterflat'''
    '''this is only used for the photometric imgs processing'''

    def Flat(self):
        name = re.split('/', self.outdir)[-2:]
        name = str(name[0])+'_'+str(name[1])+'.fits'
        print('name : ', name)

        calibratedIn = os.path.join(self.outdir, 'calibrated.fits')
        calibratedOut = calibratedIn#os.path.join(self.outdir, name)
        filt = self.hdr['FILTERS']
        logging.info('filter '+filt)
        if (filt.find('CBPG') >= 0):
            master_flat = None
            print('currently not using CBPG masterflat because gradient is much higher than in actual sky images')
            #master_flat = os.path.join(os.environ['MASTER_PREFIX'],\
            #                       'master_g_2017_10.fits')
        elif (filt.find('r') >= 0):
            master_flat = os.path.join(os.environ['MASTER_PREFIX'],
                                       'master_r_2017_10.fits')
        elif (filt.find('i') >= 0):
            master_flat = os.path.join(os.environ['MASTER_PREFIX'],
                                       'master_i_2017_10.fits')

        if master_flat:
            cmd = 'compute.py %s / %s %s %s'\
                  % (calibratedIn, '1.', master_flat, calibratedOut)
            logging.info(cmd)
            os.system(cmd)
        else:
            print('no masterflat found for this image')
        return name

    '''Use output catalog of Sextractor to find the brighter rounder object'''

    def FindObject(self):
        objects = tb.ParseSexCat(os.path.join(self.outdir, 'se.list'))
        objects = objects[(np.sqrt(objects.field('X2_IMAGE')) / np.sqrt(objects.field('Y2_IMAGE')) < 10)
                          & (np.sqrt(objects.field('X2_IMAGE')) / np.sqrt(objects.field('Y2_IMAGE')) > .1)]
        objects = objects[objects.field('A_IMAGE') < 100]
        objects = objects[objects.field('FLUX_BEST') == max(objects.field('FLUX_BEST'))]
        logging.info('x_object, y_object = ' + str(objects.field('X_IMAGE')[0]) +
                     ' ' + str(objects.field('Y_IMAGE')[0]))
        logging.debug('Object selected :' + str(objects))
        return objects.field('X_IMAGE')[0], objects.field('Y_IMAGE')[0]

    ''' run sectractor on image to extract a background map and a se.list'''

    def RunSEx(self, image):
        default = os.path.join(os.environ['SEX_PREFIX'], 'default.sex')
        cmd = "sex -c %s %s  -CATALOG_NAME=%s"%(default, image,
                                                os.path.join(self.outdir, 'se.list'))
        print(cmd)
        os.system(cmd)
        os.system('mv segmentation.fits %s' % (self.outdir))
        return


'''--------------------------------------------'''
''' A class to extract spectrum from fits image'''
'''--------------------------------------------'''


class Process(object):
    def __init__(self, Image,
                 outdir='./',
                 orientation='x',
                 plot=False,
                 dispersion=False,
                 calibrate=False,
                 Map=False,
                 method=1,
                 isr=False,
                 cosmic=False,
                 aperture=150,
                 object_name=None, **kwargs):

        self.method = method
        self.aperture = aperture
        self.calibrate = calibrate
        self.dispersion = dispersion
        self.orientation = orientation
        self.plot = plot
        self.img = Image
        self.Map = Map
        self.isr = isr
        self.cosmic = cosmic
        self.instrument = instru.telinst(self.img)
        self.object_name = object_name
        if self.object_name is None:
            self.object_name = str(self.instrument.header.get('OBJECT'))
        logging.info('Object name :' + str(self.object_name))
        self.outdir = outdir

    def SetFootprint(self, x_object, y_object,
                     order='m+1',
                     offset=100, **kwargs):
        extend = 700
        half_width = int(self.aperture/2.)
        translate_x = offset
        translate_y = offset
        disp_res = str(self.instrument.header.get('FILTERS'))
        logging.info("FILTERS : " + disp_res)
        if (disp_res.find('RONCHI400') >= 0):
            extend = 1000
        if(order == 'm-1'):
            translate_x = - extend - offset
            translate_y = - extend - offset

        if (self.orientation == 'x'):
            x_start = x_object + translate_x
            x_end = x_start+extend
            #''' footprint'''
            y_start = y_object - half_width
            y_end = y_object+half_width
            #''' footprint'''
        elif(self.orientation == 'y'):
            x_start = x_object - half_width
            x_end = x_object+half_width
            #''' footprint'''
            y_start = y_object + transalte_y
            y_end = y_start+extend
            #''' footprint'''
        logging.debug('Is ' + str(x_end)+' is smaller than '+str(self.instrument.lenX) + ' ? -> choose min')
        x_end = min(x_end, self.instrument.lenX)           #''' checkbound'''
        y_end = min(y_end, self.instrument.lenY)           #''' checkbound'''
        y_start = max(y_start, 0)                            #''' checkbound'''
        x_start = max(x_start, 0)                            #''' checkbound'''

        if(x_end < x_start):
            sys.exit('Looking for a footprint out of imager bounds, abort.\n Possibly object identification went wrong ?')
        logging.debug('(x_start, x_end) = (' + str(x_start)+', '+str(x_end)+')')
        logging.debug('(y_start, y_end) = (' + str(y_start)+', '+str(y_end)+')')

        return x_start, x_end, y_start, y_end

    '''handling the extraction of the spectrum'''

    def SpectrumExtraction(self, x_object, y_object,
                           seeing=None,
                           sx=None, sy=None,
                           orders=['m+1'],
                           offset=100, **kwargs):

        for order in orders:
            logging.info('order : ' + order)
            x_start, x_end, y_start, y_end = self.SetFootprint(x_object, y_object,
                                                               offset=offset,
                                                               order=order)

            mask = os.path.join(self.outdir, 'segmentation.fits')
            spec = spectrum.spectrum(self.object_name, x_object, y_object,
                                     os.path.join(self.outdir, 'calibrated.fits'),
                                     x_start, x_end, y_start, y_end,
                                     mask=mask,
                                     dispersion=self.dispersion,
                                     calibrate=self.calibrate,
                                     orientation=self.orientation,
                                     plot=self.plot,
                                     out_dir=self.outdir,
                                     Map=self.Map,
                                     seeing=seeing,
                                     sx=sx,
                                     sy=sy,
                                     order=order,
                                     method=self.method,
                                     isr=self.isr,
                                     cosmic=self.cosmic,
                                     offset=offset)

            '''handle the extraction of the spectrum and then compress the image'''
            direction = 'y'   # Profile direction of background projection
            spec.TargetSpectrum(direction)

        calibrated = os.path.join(self.outdir, 'calibrated.fits')
        logging.info('gzip -f '+calibrated+' and segmentation.fits')
        os.system('gzip -f %s'%(calibrated))
        os.system('gzip -f %s'%(os.path.join(self.outdir, 'segmentation.fits')))

        return spec.pixel, spec.wavelength, spec.aper_flux, spec.calib_profile


# ===========================================================
class filters(Exception):
    pass


# ===========================================================
'''A class to implement various filters on the pixels'''


class filters(object):
    def __init__(self, fitsimage=None, plot=False, *args, **kwargs):
        if fitsimage:
            (self.filepath, self.fitsimage) = os.path.split(fitsimage)
        self.plot = plot

    def Cosmics2(self, data, seeing, mean, sigma):
        Iter = 0
        count = 1000000
        lenX = len(data)
        lenY = len(data[0])
        CosmicImage = np.zeros([lenX, lenY])
        total = 0
        logging.info('Looking for cosmics with (Seeing, BkgdMean, BkgdSigma) = '
                     + str(seeing)+' ' + str(mean) + ' ' + str(sigma))
        while ((count) and (Iter < 5)):
            print(" Iter ", Iter+1)
            count = self.LaplacianFilter(sigma, mean, seeing, data, CosmicImage)
            print(" Number of cosmic found ", count)
            total += count
            Iter += 1
        return CosmicImage

    def Cosmics(self, seeing, frame=None, **kwargs):
        data = (pf.open(os.path.join(self.filepath, self.fitsimage)))[0].data
        if frame is not None:
            xi = frame[0]
            xf = frame[1]
            yi = frame[2]
            yf = frame[3]
            data = data[yi:yf, xi:xf]
            hdu = pf.PrimaryHDU(data)
            hdu.writeto(os.path.join(self.filepath, 'frame_'+self.fitsimage), overwrite=True)
        Iter = 0
        count = 1000000
        lenX = len(data)
        lenY = len(data[0])
        print('x_size, y_size : ', lenX, lenY)
        mean, sigma = tb.BkgdStatClip(data, self.plot)
        CosmicImage = np.zeros([lenX, lenY])
        total = 0
        while ((count) and (Iter < 5)):
            print(" Iter ", Iter+1)
            count = self.LaplacianFilter(sigma, mean, seeing, data, CosmicImage)
            print(" Number of cosmic found ", count)
            total += count
            Iter += 1

        ''' Write mask if cosmics found on frame'''
        if total:
            hdu = pf.PrimaryHDU(CosmicImage)
            hdu.writeto(os.path.join(self.filepath, 'cosmics_'+self.fitsimage), overwrite=True)
        return


#//Laplacian filter
#/*!Cuts (based on the article -> astro-ph/0108003):
#  -cut_lap : the laplacian operator increases the noise by a factor of
#  "sqrt(1.25)"
#
#  -cut_f : 2*sigma(med), where sigma(med) is the variance of the
#  sky's median calculated in a box (3*3), here.
#  (sigma(med) = sigma(sky)*1.22/sqrt(n); n = size of the box)
#
#  -cut_lf : calculated from the article.
#  Factor 2.35 -> to have the seeing in arc sec */


    def LaplacianFilter(self, Sigma, Mean, seeing, data, CosmicImage):
        xmax = len(data)-1
        ymax = len(data[0])-1
        l = 0.0
        f = 0.0
        med = 0.0
        cut = 4 * Sigma
        cut_lap = cut * np.sqrt(1.25)
        cut_f = 2 * (Sigma*1.22/3)
        cut_lf = 2./(seeing*2.35-1)
        count = 0
        for j in range(1, ymax):
            for i in range(1, xmax):
                #Calculation of the laplacian and the median only for pixels > 3 sigma
                if (data[i, j] > cut+Mean):
                    l = data[i, j] - 0.25*(data[i-1, j] + data[i+1, j] +
                                           data[i, j-1] + data[i, j+1])
                    med = np.median(data[i-1:i+2, j-1:j+2])
                    f = med - Mean  #f is invariant by addition of a constant
                    #Construction of a cosmic image
                    if((l > cut_lap) and ((f < cut_f) or ((l/f) > cut_lf))):
                        CosmicImage[i, j] = 1
                        #print "replacing ", i, j,  data[i,j], ' with '
                        data[i, j] = med
                        #print med
                        count += 1
                # Image is set to 0 by default
        #CosmicImage.Simplify(0.5)
        return count
