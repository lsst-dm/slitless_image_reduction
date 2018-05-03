#!/usr/bin/env python 


import os, sys
import numpy as np
import logging




# ===========================================================
class spectrum(Exception):
    pass
# ===========================================================


''' A class to extract a spectrum from a fits image'''
''' Warning : x and y from usual ds9 are exchanged due to python convention !!!'''
class dispersion(object):
    def __init__(self, run,  filters, offset = 100, **kwargs):
        self.offset        = offset
        self.run           = run
        self.filters       = filters
        self.pix2wgth_coef = self.Pixel2WavelengthCoefficients()


    ''' Find pixel2wavelength transfo, write coefficients'''
    def Pixel2WavelengthCoefficients(self):
        obs_lines = self.ObjectObservedLines()
        ref_lines = self.ReferenceLines()
        if((obs_lines is None) or (ref_lines is None)):
            obs_lines = [2,3]
            ref_lines = [2,3]
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
            if (self.filters.find('RONCHI200') >= 0):
                obs_lines = [117, 190]  #guessing half...      
        else :
            logging.info('No lines observed for the run')
            return None
        return obs_lines

    ''' spectral features (nm)'''
    def ReferenceLines(self):
        if(self.run==1):
            ref_lines = [486, 761]
        elif(self.run==2):
            ref_lines = [685, 761]
        elif(self.run==3):
            ref_lines = [486.2, 656.5, 761]
        elif(self.run==5):
            ref_lines = [486.2, 656.5, 761]
        elif(self.run==7):
            ref_lines = [486.2, 761]
        else :
            logging.info('No ReferenceLines for the run')
            return None
        return ref_lines


     
    '''Currently linear transfo'''
    def Wavelength2Pixel(self, wgth):
        wgth = np.asarray(wgth, dtype=np.float64)
        coef = self.pix2wgth_coef
        return (wgth-coef[1]) / coef[0]

    def Pixel2Wavelength(self, pixel):
        pixel = np.asarray(pixel, dtype=np.float64)
        coef = self.pix2wgth_coef
        return np.array(coef[1] + coef[0] * pixel)




