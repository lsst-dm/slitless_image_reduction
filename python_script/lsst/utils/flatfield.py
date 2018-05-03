#!/usr/bin/env python 

import os, sys
import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate as interp
from scipy.interpolate import griddata
import toolbox as tb
import logging
from astropy.stats import sigma_clip
import extraction as ex
import telinst as instru
import astropy.io.fits as pf


'''this class must be called after spectrum class as it extensively '''
'''uses Pixel2WavelengthCoefficients()'''
class flatfield(object):
    def __init__(self, spectrum, footprint, head, dispers, **kwargs):
        self.footprint  = footprint
        self.spectrum   = spectrum
        self.dispersion = dispers
        self.head       = head
        self.FlatList()
   

    '''if able to buid a map of pixel defect from flats-> divide footprint by it'''
    def divideByMap(self, footprint):
        defect_map = self.BuildDefectMap()
        hdu = pf.PrimaryHDU(defect_map)
        name = str("defect_map"+self.spectrum.order+".fits")
        hdu.writeto(os.path.join(self.spectrum.out_dir,name), overwrite=True)
        logging.info("dividing footprint by map of defects "+ self.masterflat_list)
 
        footprint = footprint/defect_map
        return footprint
        
        
    def FlatList(self):
        self.monochro_flat_path\
            = os.path.join(os.environ['PROCESS_LSST_AUX_TEL'],'monochromatic_flatfield/')
        #file =  'synthetic_flats_RONCHI.list'
        file =  'synthetic_flats.list'
        self.monochro_flat_list = os.path.join(self.monochro_flat_path, file)
        self.masterflat_path = os.environ['MASTER_PREFIX']+'/'
        if((self.spectrum.run == 1) or (self.spectrum.run == 2)):
            self.masterflat_list =  os.path.join(self.masterflat_path, 'masterflat.list')
        elif(self.spectrum.run == 3):
            self.masterflat_list =  os.path.join(self.masterflat_path, 'masterflatNov2016.list')
        elif(self.spectrum.run == 5):
            self.masterflat_list =  os.path.join(self.masterflat_path, 'masterflat_run5.list')
        elif(self.spectrum.run == 7):
            self.masterflat_list =  os.path.join(self.masterflat_path, 'masterflatOct2017.list')
        else:
            self.masterflat_list =  None 
        return 

        
        
    def Calibrate(self, footprint, head):
        if self.monochro_flat_list :
            logging.info("A synthetic flat is built")
            synthetic_flat = self.BuildSyntheticFlat()
        else :
            sys.exit('Should propose other pixel calibration\
            than from monochromatic flatfields')

        fig=pl.figure()
        footprint = footprint / synthetic_flat # out in Nb of photons
        #pl.imshow(footprint, cmap='hot', aspect='auto')
        #pl.show()
        X, calib_flat = self.project(synthetic_flat)          
        '''Extract after calibration'''
        get = ex.extract(self.spectrum, footprint)
        X, calib_profile  = get.ExtractSpectrum()
        '''Dump spectrum'''
        logging.info("Writing calibrated spectrum using synthetic flat into calibspectrum.list")
        self.spectrum.DumpSpectrum(head,
                          ('pixel', 'w','rawflux', 'califlux'),
                          (self.spectrum.pixel, self.spectrum.wavelength,
                           np.array(self.spectrum.aper_flux), np.array(calib_profile)),
                          str("calibspectrum"+self.spectrum.order+".list") )           
   
        if self.spectrum.plot :
            self.ControlPlot1(X, calib_profile, calib_flat, self.spectrum.aper_flux)
       

    def BuildSyntheticFlat(self):
        ''' monitoring photodiode, R is in (mA/W)''' 
        wg, responsivity, qe =\
        tb.readfile(self.spectrum.monitoring_photodiode_file ,\
                    ['wavelength', 'responsivity','QE'])

        xsize         = self.spectrum.x_end - self.spectrum.x_start
        ysize         = self.spectrum.y_end - self.spectrum.y_start
        syntheticflat = np.zeros([ysize, xsize])
        pointsX =[] ; pointsY = [] ; values = []    
        fits, wgth   = tb.readtxt(self.monochro_flat_list, ['fits', 'wavelength'])
        for img, w in zip(fits, wgth):
            pixel      = self.dispersion.Wavelength2Pixel(w)
            # ! pixel is relative to the object-> translate with respect to frame
            pixel      = self.dispersion.RefPos2FramePos(pixel)
            flat       = instru.telinst(self.monochro_flat_path + img)
            head       = flat.header  
            iphtot     = head['IPHTOT'] 
            siphtot    = head['SIPHTOT']
            image_data = tb.Flip((flat.Image(self.monochro_flat_path + img))\
                              [self.spectrum.y_start:self.spectrum.y_end, self.spectrum.x_start:self.spectrum.x_end], flip=self.spectrum.flip)   
            which_amp  = flat.IsInAmp(self.spectrum.y_start, self.spectrum.x_start, Trimmed=True)
            low        = head['LOW'+str(which_amp)]  
            high       = head['HIGH'+str(which_amp)]
            '''remove 5s outliers'''
            if(self.spectrum.orientation=='x'):
                signal     = tb.Flip((image_data[:, int(pixel)])[(image_data[:, int(pixel)] > low)\
                                            & (image_data[:, int(pixel)] < high)],flip=self.spectrum.flip)
            if(self.spectrum.orientation=='y'):
                signal     = tb.Flip((image_data[int(pixel),:])[(image_data[int(pixel),:] > low)\
                                            & (image_data[int(pixel),:] < high)],flip=self.spectrum.flip)
     
            respinterp = interp.griddata(wg, responsivity, w)     
            signal     =   np.mean(signal) / (iphtot/(respinterp))
            '''Convert flux into nb of e_ because CCD counts e-, then to QE''' 
            logging.info('calib : '+str(w)+str(signal))
            
            logging.debug('Normalizing flat '+ img+' by photodiode QE')
            if(self.spectrum.plot):
                if(self.spectrum.orientation=='x'):
                    syntheticflat[:, int(pixel)] = signal 
                if(self.spectrum.orientation=='y'):
                    syntheticflat[int(pixel), :] = signal 

            for i in range(int(ysize)):
                pointsX.append(int(pixel))
                pointsY.append(i)
                values.append(signal)

        values = values / max(values)
        ''' Build the synthetic flat'''
        points         = np.vstack((np.array(pointsY), np.array(pointsX))).T
        values         = (np.array(values)).T
        grid_x, grid_y = np.mgrid[0:ysize, 0:xsize]
        interp_map     = griddata(points, values, (grid_x, grid_y), method='linear')

        '''Control plots'''
        if(self.spectrum.plot):
            self.ControlPlot2(values, syntheticflat, interp_map)
            
        return interp_map

   



    

    '''To flatfield the pixel defects '''
    def BuildDefectMap(self):
        xsize         = int(np.round(self.spectrum.x_end - self.spectrum.x_start,3))
        ysize         = int(np.round(self.spectrum.y_end - self.spectrum.y_start,3))
        defect_map    = np.ones([ysize, xsize])
   
        if self.spectrum.orientation == 'y' :
            defect_map = defect_map.T
            xsize      = ysize
        '''Bad S/N for mono-> switch to broadband flat'''
        #fits, wgth  = tb.readtxt(self.monochro_flat_list, ['fits', 'wavelength'])
        #path = self.monochro_flat_path
        fits, wgth  = tb.readtxt(self.masterflat_list, ['fits', 'wavelength'])
        out         = zip(*sorted(zip(wgth, fits)))
        length      = len(out[0])
        if length == 1 :
            logging.info("using flat: "+ self.masterflat_path + out[1][0])
            ref1  = instru.telinst(self.masterflat_path + out[1][0])
            defect_map = tb.Flip((ref1.Image(self.masterflat_path + out[1][0]))\
                              [self.spectrum.y_start : self.spectrum.y_end, self.spectrum.x_start : self.spectrum.x_end], flip=self.spectrum.flip) 
        elif length > 1 :
            defect_map = self.SyntheticMap(defect_map, length, out, xsize)
            if ((self.spectrum.orientation == 'y') and(self.defect_map )) :
                defect_map = defect_map.T # To come back to the footprint
        return defect_map
    

    
    '''Interpolating between various flat to build a map of defect evolving with wght'''
    def SyntheticMap(self, defect_map, length, out, xsize):
        path   = self.masterflat_path
        logging.info("Building a synthetic flat from "+ str(path) + str(out[1][:]))
        for i in range(length-1):
            ref1  = instru.telinst(path + out[1][i])
            img1  = tb.Flip((ref1.Image(path + out[1][i]))\
                         [self.spectrum.y_start : self.spectrum.y_end, self.spectrum.x_start : self.spectrum.x_end],flip=self.spectrum.flip)
            ref2  = instru.telinst(path + out[1][i+1])
            img2  = tb.Flip((ref2.Image(path + out[1][i+1]))\
                         [self.spectrum.y_start : self.spectrum.y_end, self.spectrum.x_start : self.spectrum.x_end], flip=self.spectrum.flip)


            if self.spectrum.orientation == 'y': # to process in the direction needed
                img1  = img1.T
                img2  = img2.T

            wgth1 = out[0][i]
            wgth2 = out[0][i+1]
            # ! pixel is relative to the object-> translate with respect to frame
            pixel1 = int(self.spectrum.RefPos2FramePos(self.dispersion.Wavelength2Pixel(wgth1)))
            pixel2 = int(self.spectrum.RefPos2FramePos(self.dispersion.Wavelength2Pixel(wgth2)))

            if(pixel1 >= self.spectrum.x_end - self.spectrum.x_start):             
                logging.info('Defect map out of bonds : the footprint found is likely mostly outside of frame ->  divide by ones')
                return defect_map
            for position in range(pixel1, pixel2):                
                norm   = np.median(img1[:, position])
                row1   = img1[:, position]/norm
                norm   = np.median(img2[:, position])
                row2   = img2[:, position]/norm
                coef   = (float(position)-float(pixel1)) / (float(pixel2)-float(pixel1))
                defect_map[:, position]  = row1 * (1-coef)\
                                          + row2 * coef
                if(i==0): # Extending map to the entire footprint
                    for position in range(0, pixel1):
                        norm = np.median(img1[:, position])
                        row  = img1[:, position]/norm
                        defect_map[:, position]  = row
                if(i+1==(length-1)):
                    for position in range(pixel2, xsize):
                        norm = np.median(img2[:, position])
                        row  = img2[:, position]/norm
                        defect_map[:, position]  = row
        return defect_map



    '''  Extracting spectrum from array     ''' 
    def project(self, footprint):
        pixel = [] ; aper_flux = []
        if self.spectrum.orientation == 'y':
            footprint  = footprint.T
        extend    = len(footprint)
        pixels    = np.arange(0, extend)
        end       = len(footprint[0])
        for i in range(end):
            projected_footprint = footprint[: , i]
            pixel.append(i)
            '''Aperture flux'''
            Sum  = np.sum(projected_footprint)
            aper_flux.append(Sum)      
        return pixel, aper_flux
    




    def ControlPlot1(self, pixel, calib_profile, calib_flat, raw_profile):
        # ! Pixel is relative to frarme->translate with respect to Ref
        # ! so that transfo applies.
        wX  = self.dispersion.FramePos2RefPos(pixel)
        wX  = self.dispersion.Pixel2Wavelength(wX)                         
        fig = pl.figure(3)
        pl.plot(wX,calib_profile, color='magenta', label='normalized spectrum' )
        pl.legend()
        fig.savefig(os.path.join(self.spectrum.out_dir,"normalized_spectrum.pdf"))
        fig = pl.figure(4)
        pl.plot(wX,calib_flat, color='magenta', label='synthetic flat profile')
        pl.legend()
        fig.savefig(os.path.join(self.spectrum.out_dir,"synth_flat_profile.pdf"))
        fig = pl.figure(5)
        pl.plot(wX, raw_profile, color='magenta', label= 'raw spectrum')
        pl.legend()
        fig.savefig(os.path.join(self.spectrum.out_dir,"raw_spectrum.pdf"))
        return


    def ControlPlot2(self, values, syntheticflat, interp_map):
        fig=pl.figure(7,figsize=(8, 4)) 
        pl.subplot(1,2,1)
        pl.imshow(syntheticflat/max(values), cmap='hot', aspect='auto')
        pl.ylabel('pixel')
        pl.xlabel('pixel')
        pl.subplot(1,2,2)
        pl.imshow(interp_map/max(values), cmap='hot', aspect='auto')
        pl.xlabel('pixel')
        pl.suptitle('Synthetic flat (450:950)nm')
        pl.colorbar()
        fig.savefig(os.path.join(self.spectrum.out_dir,'synthflat.pdf'))
        return

