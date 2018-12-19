#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import next
from builtins import str
from builtins import zip
from builtins import range
from builtins import object
import telinst as instru
import matplotlib.cm as cm
import os
import sys
import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as pl
import croaks
import scipy.interpolate as interp
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy import integrate
import operator
from croaks import NTuple
import toolbox as tb
import astropy.time
import dateutil.parser
import reduceimage as ri
import extraction as ex
import logging
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def Flip(frame, flip=False, **kwargs):
    if flip == True:
        frame = frame[:, ::-1]
    return frame


# Define model function to be used to fit to the data above:
def Gauss1D(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def MoffatFit(pixels, projected_footprint, A, mu, sigma):
    g_init = models.Moffat1D(amplitude=A, x_0=mu, gamma=sigma)
    fit_g = fitting.LevMarLSQFitter()
    psf = fit_g(g_init, pixels, projected_footprint)
    start = psf.x_0 - 5 * psf.gamma
    end = psf.x_0 + 5 * psf.gamma
    integral = (integrate.quad(lambda pixels: psf(pixels), start, end))[0]
    ''' begin Control plot'''
    #pl.plot(pixels, psf(pixels), label='Moffat')
    #pl.yscale('log')
    #pl.ylim(1., 1E6)
    #pl.plot(projected_footprint)
    ''' End control plot'''
    return


def VoigtFit(pixels, profile, A, mu, sigma, weights=[], **kwargs):
    #pixels = pixels[20:60] est a smaller footprint. bof. definitley need a weight map.
    #profile = profile[20:60]
    ### I will return the gaussian profile as a weight map
    #weights=1/profile**2
    weights = np.ones(len(pixels))
    g_init = models.Voigt1D(amplitude_L=A, x_0=mu, fwhm_G=sigma*2., fwhm_L=sigma*3)
    #g_init.x_0.min = 0.
    #g_init.x_0.max = 80.
    #g_init.fwhm_L.min = 0.1
    #g_init.fwhm_G.min = 0.1
    #g_init.fwhm_L.max = 20.
    #g_init.fwhm_G.max = 20.
    fit_g = fitting.LevMarLSQFitter()
    psf = fit_g(g_init, pixels, profile, weights=weights)
    start = psf.x_0 - 5 * (psf.fwhm_L+psf.fwhm_G)/2
    end = psf.x_0 + 5 * (psf.fwhm_L+psf.fwhm_G)/2
    integral = (integrate.quad(lambda pixels: psf(pixels), start, end))[0]
    ''' begin Control plot'''
    #pl.plot(pixels, psf(pixels), label='Voigt')
    #pl.plot(pixels, weights, label='weights')
    #pl.yscale('log')
    #pl.ylim(1., 1E6)
    #pl.plot(profile)
    #pl.legend()
    #pl.show()
    ''' End control plot'''
    return integral, psf.fwhm_L.value, psf.fwhm_G.value


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
                 order='m+1',
                 method=1,
                 offset=100, **kwargs):

        self.aperture = True #So far, only this choice implemented
        self.calibrate = calibrate
        self.orientation = orientation
        self.dispersion = dispersion
        self.out_dir = out_dir
        self.plot = plot
        self.Map = Map
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
        self.FlatList()
        self.pix2wgth_coef = [0]
        if self.dispersion is True:
            self.pix2wgth_coef = self.Pixel2WavelengthCoefficients()
        self.aper_flux = None
        self.calib_profile = None
        self.pixel = None
        self.wavelength = None
        if (self.order == 'm+1'):
            self.flip = False
        elif(self.order == 'm-1'):
            self.flip = True

    def FlatList(self):
        self.monochro_flat_path\
            = os.path.join(os.environ['PROCESS_LSST_AUX_TEL'], 'monochromatic_flatfield/')
        #file =  'synthetic_flats_RONCHI.list'
        file = 'synthetic_flats.list'
        self.monochro_flat_list = os.path.join(self.monochro_flat_path, file)
        self.masterflat_path = os.environ['MASTER_PREFIX']+'/'
        if((self.run == 1) or (self.run == 2)):
            self.masterflat_list = os.path.join(self.masterflat_path, 'masterflat.list')
        if(self.run == 3):
            self.masterflat_list = os.path.join(self.masterflat_path, 'masterflatNov2016.list')
        if(self.run == 5):
            self.masterflat_list = os.path.join(self.masterflat_path, 'masterflat_run5.list')
        return

    ''' This is used to model the periodic pattern\
    imprint by the electronic on the background'''

    def ExtractProfile(self, footprint, mask, direction):
        Y = []
        if(direction == "y"):
            footprint = footprint.T
            mask = mask.T
        rows, cols = len(footprint), len(footprint[0])
        for i in range(0, rows):
            keep = footprint[i, :]
            remo = mask[i, :]
            keep = keep[remo == 0]
            if len(keep) is 0:
                flux = Y[-1]
            else:
                flux = np.median(keep)
            Y.append(flux)

        rows = np.arange(0, rows)
        file = str("background"+self.order+".list")
        tb.DumpTuple(('rows', 'level'), (rows, Y), os.path.join(self.out_dir, file))
        return Y

    '''Catch surrounding frame of footprint and masks from mask'''
    ''' Also return footprint's mask                           '''

    def Mask(self, inst, direction):
        which_amp = inst.IsInAmp(self.y_start, self.x_start, Trimmed=True)
        amp_frame = inst.TelAmpLimits(which_amp, Trimmed=True)
        masking = instru.telinst(self.mask, verbose='')
        m = masking.Image(self.mask)
        footprint_mask = Flip(m[self.y_start: self.y_end,
                                self.x_start: self.x_end], flip=self.flip)  # Select the footprint mask
        '''A surrounding frame to extract the periodic pattern in the background'''
        y_mask_start = max(self.y_start-100, amp_frame[0]) # x<->y
        y_mask_end = min(self.y_end+100, amp_frame[1]) # x<->y
        m = Flip(m[y_mask_start:y_mask_end, self.x_start:self.x_end], flip=self.flip)
        f = Flip((inst.Image(self.image))[y_mask_start:y_mask_end, self.x_start:self.x_end],
                 flip=self.flip)

        return self.ExtractProfile(f, m, direction), footprint_mask

    '''The actual extraction of the spectrum'''

    def TargetSpectrum(self, Direction='y', **kwargs):
        logging.info('Opening of image : ' + str(self.image))
        inst = instru.telinst(self.image, verbose='')
        image_data = inst.Image(self.image)
        head = inst.header
        name = head['OBJECT'] + head['RECID']
        saturation = np.max(image_data)*0.8                # get saturation
        footprint = Flip(image_data[self.y_start: self.y_end,
                                    self.x_start: self.x_end], flip=self.flip)  # Select the footprint

        '''Subtract profile'''
        profile = []
        if(self.mask):
            profile, footprint_mask = self.Mask(inst, Direction)
            if (len(profile) != 0): # subtract periodic background
                footprint = self.SubtractProfile(footprint, profile)

        '''Determine and write map of defects   '''
        '''if True, Correct for defective pixels'''
        if ((self.Map is True) and (self.dispersion is True)):
            defect_map = self.BuildDefectMap()
            hdu = pf.PrimaryHDU(defect_map)
            name = str("defect_map"+self.order+".fits")
            hdu.writeto(os.path.join(self.out_dir, name), overwrite=True)
            footprint = footprint/defect_map

        ''' Write footprint before removing cosmics in case of other than aperture phot'''
        if not self.aperture:
            hdu = pf.PrimaryHDU(footprint)
            name = str("footprint"+self.order+".fits")
            hdu.writeto(os.path.join(self.out_dir, name), overwrite=True)

        '''Extract a map of cosmic in footprint'''
        if(self.mask):
            mean, sigma = tb.BkgdStat(footprint, footprint_mask)
            cosmics = ri.filters(plot=self.plot)
            cosmicimage = cosmics.Cosmics2(footprint, self.seeing, mean, sigma)
            hdu = pf.PrimaryHDU(cosmicimage)
            name = str("cosmics"+self.order+".fits")
            hdu.writeto(os.path.join(self.out_dir, name), overwrite=True)

        ''' Write footprint after removing cosmics in case of aperture phot'''
        if self.aperture:
            hdu = pf.PrimaryHDU(footprint)
            name = str("footprint"+self.order+".fits")
            hdu.writeto(os.path.join(self.out_dir, name), overwrite=True)

        ''' Extract raw profile'''
        ''' Testing two versions:'''
        '''1st is constant aperture'''
        '''second used gaussian fit to determine center and width'''
        logging.info("Spectrum extraction method : "+str(self.method))

        get = ex.extract(self, footprint)

        if (self.method == 1):
            #self.pixel, self.aper_flux = self.ExtractSpectrum(footprint)
            self.pixel, self.aper_flux = get.ExtractSpectrum()

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
            self.wavelength = self.Pixel2Wavelength(self.pixel)
            '''normalize using synthetic flat'''
            if self.calibrate is True:
                self.Calibrate(footprint, head)
        else:
            self.wavelength = np.zeros(len(self.pixel))

        '''Dump raw spectrum'''
        self.DumpSpectrum(head,
                          ('pixel', 'w', 'aper_flux',
                           'psf_gauss_flux', 'psf_gauss_sigma', 'psf_gauss_mu',
                           'psf_voigt_flux', 'psf_voigt_fwhmL', 'psf_voigt_fwhmG'),
                          (self.pixel, self.wavelength, np.array(self.aper_flux),
                           np.array(get.psf_gauss_flux),
                           np.array(get.psf_gauss_sigma),
                           np.array(get.psf_gauss_mu),
                           np.array(get.psf_voigt_flux),
                           np.array(get.psf_voigt_fwhmL),
                           np.array(get.psf_voigt_fwhmG)),
                          str("rawspectrum"+self.order+".list"))

        return

    def Calibrate(self, footprint, head):
        if self.monochro_flat_list:
            logging.info("A synthetic flat is built")
            synthetic_flat = self.BuildSyntheticFlat()
        else:
            sys.exit('Should propose other pixel calibration\
            than from monochromatic flatfields')

        fig = pl.figure()
        footprint = footprint / synthetic_flat # out in Nb of photons
        #pl.imshow(footprint, cmap='hot', aspect='auto')
        #pl.show()
        X, calib_flat = self.project(synthetic_flat)
        '''Extract after calibration'''
        get = ex.extract(self, footprint)
        X, self.calib_profile = get.ExtractSpectrum()
        #X, self.calib_profile = self.ExtractSpectrum(footprint)
        '''Dump spectrum'''
        self.DumpSpectrum(head,
                          ('pixel', 'w', 'rawflux', 'califlux'),
                          (self.pixel, self.wavelength,
                           np.array(self.aper_flux), np.array(self.calib_profile)),
                          str("calibspectrum"+self.order+".list"))

        if self.plot:
            self.ControlPlot1(X, self.calib_profile, calib_flat, self.aper_flux)

    def DumpSpectrum(self, head, names, list, file):
        list = list(zip(*list))
        info = np.rec.fromrecords([i for i in list], names=names)

        info = info.view(NTuple)
        if (self.filters.find('RONCHI400') >= 0):
            info.keys['RONCHI'] = 400
        else:
            info.keys['RONCHI'] = 200
        info.keys['IMG_NAME'] = str(head.get('IMG_NAME')+' ' + self.out_dir)
        info.keys['AIRMASS'] = head.get('AIRMASS')
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
        info.keys['PIX2WGTH'] = self.pix2wgth_coef[0]
        info.keys['FLATFIELDED'] = self.calibrate
        info.keys['DISPERSION'] = self.dispersion
        info.keys['ORIENTATION'] = self.orientation
        info.keys['SEEING'] = self.seeing
        info.keys['SX'] = self.sx
        info.keys['SY'] = self.sy
        info.keys['OBJECT'] = self.object_name
        info.keys['INFO1'] = "pixel # start from direct object"
        info.totxt(os.path.join(self.out_dir, file))
        return

    '''subtract a longitudinal profile to the footprint'''

    def SubtractProfile(self, footprint, Profile):
        for i in range(len(Profile)):
            footprint[:, i] = footprint[:, i] - Profile[i]
        return footprint

    '''Interpolating between various flat to build a map of defect evolving with wght'''

    def SyntheticMap(self, defect_map, length, out, xsize):
        path = self.masterflat_path
        for i in range(length-1):
            ref1 = instru.telinst(path + out[1][i])
            img1 = Flip((ref1.Image(path + out[1][i]))
                        [self.y_start: self.y_end, self.x_start: self.x_end], flip=self.flip)
            ref2 = instru.telinst(path + out[1][i+1])
            img2 = Flip((ref2.Image(path + out[1][i+1]))
                        [self.y_start: self.y_end, self.x_start: self.x_end], flip=self.flip)

            if self.orientation == 'y': # to process in the direction needed
                img1 = img1.T
                img2 = img2.T

            wgth1 = out[0][i]
            wgth2 = out[0][i+1]
            # ! pixel is relative to the object-> translate with respect to frame
            pixel1 = int(self.RefPos2FramePos(self.Wavelength2Pixel(wgth1)))
            pixel2 = int(self.RefPos2FramePos(self.Wavelength2Pixel(wgth2)))
            for position in range(pixel1, pixel2):
                norm = np.median(img1[:, position])
                row1 = img1[:, position]/norm
                norm = np.median(img2[:, position])
                row2 = img2[:, position]/norm
                coef = (float(position)-float(pixel1)) / (float(pixel2)-float(pixel1))
                defect_map[:, position] = row1 * (1-coef)\
                    + row2 * coef
                if(i == 0): # Extending map to the entire footprint
                    for position in range(0, pixel1):
                        norm = np.median(img1[:, position])
                        row = img1[:, position]/norm
                        defect_map[:, position] = row
                if(i+1 == (length-1)):
                    for position in range(pixel2, xsize):
                        norm = np.median(img2[:, position])
                        row = img2[:, position]/norm
                        defect_map[:, position] = row
        return defect_map

    def BuildDefectMap(self):
        xsize = int(np.round(self.x_end - self.x_start, 3))
        ysize = int(np.round(self.y_end - self.y_start, 3))
        defect_map = np.zeros([ysize, xsize])
        if self.orientation == 'y':
            defect_map = defect_map.T
            xsize = ysize
        '''Bad S/N for mono-> switch to broadband flat'''
        #fits, wgth  = tb.readtxt(self.monochro_flat_list, ['fits', 'wavelength'])
        #path = self.monochro_flat_path
        fits, wgth = tb.readtxt(self.masterflat_list, ['fits', 'wavelength'])
        out = list(zip(*sorted(zip(wgth, fits))))
        length = len(out[0])
        if length == 1:
            ref1 = instru.telinst(self.masterflat_path + out[1][0])
            defect_map = Flip((ref1.Image(self.masterflat_path + out[1][0]))
                              [self.y_start: self.y_end, self.x_start: self.x_end], flip=self.flip)
        elif length > 1:
            defect_map = self.SyntheticMap(defect_map, length, out, xsize)
            if self.orientation == 'y':
                defect_map = defect_map.T # To come back to the footprint
        return defect_map

    def BuildSyntheticFlat(self):
        ''' monitoring photodiode, R is in (mA/W)'''
        wg, responsivity, qe =\
            tb.readfile(self.monitoring_photodiode_file,
                        ['wavelength', 'responsivity', 'QE'])

        xsize = self.x_end - self.x_start
        ysize = self.y_end - self.y_start
        syntheticflat = np.zeros([ysize, xsize])
        pointsX = []
        pointsY = []
        values = []
        fits, wgth = tb.readtxt(self.monochro_flat_list, ['fits', 'wavelength'])
        for img, w in zip(fits, wgth):
            pixel = self.Wavelength2Pixel(w)
            # ! pixel is relative to the object-> translate with respect to frame
            pixel = self.RefPos2FramePos(pixel)
            flat = instru.telinst(self.monochro_flat_path + img)
            head = flat.header
            iphtot = head['IPHTOT']
            siphtot = head['SIPHTOT']
            image_data = Flip((flat.Image(self.monochro_flat_path + img))
                              [self.y_start:self.y_end, self.x_start:self.x_end], flip=self.flip)
            which_amp = flat.IsInAmp(self.y_start, self.x_start, Trimmed=True)
            low = head['LOW'+str(which_amp)]
            high = head['HIGH'+str(which_amp)]
            '''remove 5s outliers'''
            if(self.orientation == 'x'):
                signal = Flip((image_data[:, int(pixel)])[(image_data[:, int(pixel)] > low)
                                                          & (image_data[:, int(pixel)] < high)], flip=self.flip)
            if(self.orientation == 'y'):
                signal = Flip((image_data[int(pixel), :])[(image_data[int(pixel), :] > low)
                                                          & (image_data[int(pixel), :] < high)], flip=self.flip)

            respinterp = interp.griddata(wg, responsivity, w)
            signal = np.mean(signal) / (iphtot/(respinterp))
            '''Convert flux into nb of e_ because CCD counts e-, then to QE'''
            logging.info('calib : '+str(w)+str(signal))

            logging.debug('Normalizing flat ' + img+' by photodiode QE')
            if(self.plot):
                if(self.orientation == 'x'):
                    syntheticflat[:, int(pixel)] = signal
                if(self.orientation == 'y'):
                    syntheticflat[int(pixel), :] = signal

            for i in range(int(ysize)):
                pointsX.append(int(pixel))
                pointsY.append(i)
                values.append(signal)

        values = values / max(values)
        ''' Build the synthetic flat'''
        points = np.vstack((np.array(pointsY), np.array(pointsX))).T
        values = (np.array(values)).T
        grid_x, grid_y = np.mgrid[0:ysize, 0:xsize]
        interp_map = griddata(points, values, (grid_x, grid_y), method='linear')

        '''Control plots'''
        if(self.plot):
            self.ControlPlot2(values, syntheticflat, interp_map)

        return interp_map

    ''' Find pixel2wavelength transfo, write coefficients'''

    def Pixel2WavelengthCoefficients(self):
        obs_lines = self.ObjectObservedLines()
        ref_lines = self.ReferenceLines()
        if((obs_lines is None) or (ref_lines is None)):
            obs_lines = [2, 3]
            ref_lines = [2, 3]
            logging.info('No input for Pixel2WavelengthCoefficients,\
            default transformation : 1 to 1 ')
        pix2wgth_coef = np.polyfit(obs_lines, ref_lines, deg=1)

        logging.debug('Pixel -> Wavelength linear transformation coefficients : '+str(pix2wgth_coef))
        return pix2wgth_coef

    ''' Distance from the direct image (pixel)'''

    def ObjectObservedLines(self):
        if(self.run == 1):
            obs_lines = [314, 499]
        elif(self.run == 2):
            obs_lines = [324, 360]
        elif(self.run == 3):
            obs_lines = [225, 305.5, 357]
        elif(self.run == 5):
            if (self.filters.find('RONCHI400') >= 0):
                obs_lines = [450, 620, 731]
            if (self.filters.find('RONCHI200') >= 0):
                obs_lines = [223, 303, 354]
        elif(self.run == 7): #Oct2017
            if (self.filters.find('RONCHI400') >= 0):
                obs_lines = [235, 380]
        #Only two img with ronchi200, bisregard.
        else:
            logging.info('No lines observed for the run')
            return None
        return obs_lines

    ''' spectral features (nm)'''

    def ReferenceLines(self):
        if(self.run == 1):
            ref_lines = [486, 761]
        elif(self.run == 2):
            ref_lines = [685, 761]
        elif(self.run == 3):
            ref_lines = [486.2, 656.5, 761]
        elif(self.run == 5):
            ref_lines = [486.2, 656.5, 761]
        elif(self.run == 7):
            ref_lines = [486.2, 761]
        else:
            logging.info('No ReferenceLines for the run')
            return None
        return ref_lines

    ''' Position relative to direct image -> Position relative to footprint frame'''

    def RefPos2FramePos(self, pixel):
        pixel = np.asarray(pixel, dtype=np.float64)
        return pixel - self.offset

    ''' The other way around'''

    def FramePos2RefPos(self, pixel):
        pixel = np.asarray(pixel, dtype=np.float64)
        return pixel + self.offset

    '''Currently linear transfo'''

    def Wavelength2Pixel(self, wgth):
        wgth = np.asarray(wgth, dtype=np.float64)
        coef = self.pix2wgth_coef
        return (wgth-coef[1]) / coef[0]

    def Pixel2Wavelength(self, pixel):
        pixel = np.asarray(pixel, dtype=np.float64)
        coef = self.pix2wgth_coef
        return np.array(coef[1] + coef[0] * pixel)

    def ControlPlot1(self, pixel, calib_profile, calib_flat, raw_profile):
        # ! Pixel is relative to frarme->translate with respect to Ref
        # ! so that transfo applies.
        wX = self.FramePos2RefPos(pixel)
        wX = self.Pixel2Wavelength(wX)
        fig = pl.figure(3)
        pl.plot(wX, calib_profile, color='magenta', label='normalized spectrum')
        pl.legend()
        fig.savefig(os.path.join(self.out_dir, "normalized_spectrum.pdf"))
        fig = pl.figure(4)
        pl.plot(wX, calib_flat, color='magenta', label='synthetic flat profile')
        pl.legend()
        fig.savefig(os.path.join(self.out_dir, "synth_flat_profile.pdf"))
        fig = pl.figure(5)
        pl.plot(wX, raw_profile, color='magenta', label='raw spectrum')
        pl.legend()
        fig.savefig(os.path.join(self.out_dir, "raw_spectrum.pdf"))
        return

    def ControlPlot2(self, values, syntheticflat, interp_map):
        fig = pl.figure(7, figsize=(8, 4))
        pl.subplot(1, 2, 1)
        pl.imshow(syntheticflat/max(values), cmap='hot', aspect='auto')
        pl.ylabel('pixel')
        pl.xlabel('pixel')
        pl.subplot(1, 2, 2)
        pl.imshow(interp_map/max(values), cmap='hot', aspect='auto')
        pl.xlabel('pixel')
        pl.suptitle('Synthetic flat (450:950)nm')
        pl.colorbar()
        pl.legend()
        fig.savefig(os.path.join(self.out_dir, 'synthflat.pdf'))
        return

    '''  Extracting spectrum from array    mode are either [aperture, psf, both] '''

    def ExtractSpectrum2(self, footprint,
                         mode='aperture', orientation='x',
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
        if orientation == 'y':
            footprint = footprint.T

        amplitude = 1.
        extend = len(footprint)
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

        end = len(footprint[0])
        for i in range(end):
            if position is not None:
                if(min[i] < 0):
                    min[i] = 0
                projected_footprint = footprint[min[i]:max[i], i]
            else:
                projected_footprint = footprint[:, i]
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

    '''  Extracting spectrum from array     '''

    def project(self, footprint):
        pixel = []
        aper_flux = []
        if self.orientation == 'y':
            footprint = footprint.T
        extend = len(footprint)
        pixels = np.arange(0, extend)
        end = len(footprint[0])
        for i in range(end):
            projected_footprint = footprint[:, i]
            pixel.append(i)
            '''Aperture flux'''
            Sum = np.sum(projected_footprint)
            aper_flux.append(Sum)
        return pixel, aper_flux


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
        #pl.yscale('log')
    pl.xlabel('Spectrum spatial profile (pixel)')
    pl.ylabel('Amplitude (norm)')
    pl.show()
