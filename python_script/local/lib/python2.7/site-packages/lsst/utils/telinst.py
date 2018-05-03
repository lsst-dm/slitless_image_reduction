#!/usr/bin/env python 

import os, sys, re
import astropy.io.fits as pf
import numpy as np
import re
import logging
# ===========================================================
class telinst(Exception):
    pass
# ===========================================================



      
class telinst(object):
    def __init__(self, fitsimage, Level = None, **kwargs):
        logging.basicConfig(level=Level)
        
        list            = self.Acceptor(fitsimage)
        self.instrument = list[0]
        self.nccd       = list[1]
        self.namp       = list[2]
        self.lenX       = list[3]
        self.lenY       = list[4]
        self.header     = list[5]
        self.prod_path  = list[6]
        self.amp        = list[7]
        self.XYoutSize()
        
    def Acceptor(self, fitsimage):
        head = (pf.open(fitsimage))[0].header
        if ((str(head['DETECTOR']) == 'Tek2K_3') \
            and (str(head['TELID']) == 'ct36') \
            and (str(head['NAMPSYX']) == '2 2')):
            self.camera = '4_amp'          
            nccd       = 1
            namp       = 4
            lenX       = int(head['XLENGTH'])
            lenY       = int(head['YLENGTH'])
            instrument = 'ctio.9'
            prod_path  = "/Users/lpnhe/harvard/prod/"
            amp        = ('11','12','21','22')
	    return ['ctio.9', nccd, namp, lenX, lenY, head, prod_path, amp]
        elif ((str(head['DETECTOR']) == 'Tek2K_3') \
            and (str(head['TELID']) == 'ct36') \
            and (str(head['NAMPSYX']) == '1 1')):
            self.camera = 'single_amp'          
            nccd       = 1
            namp       = 1
            lenX       = 1024 # bug in int(head['XLENGTH'])
            lenY       = 1024 # bug in int(head['YLENGTH'])
            instrument = 'ctio.9'
            prod_path  = "/Users/lpnhe/harvard/prod/"
            a = str(head['AMPLIST'])
            amp        = [str(head['AMPLIST'])]
	    return ['ctio.9', nccd, namp, lenX, lenY, head, prod_path, amp]
        else :
            return NULL


    def TelAmpLimits(self, iamp=None, Trimmed=False, **kwargs):
        if (Trimmed==False):
            if (iamp=='11'):
                axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC11']))  
            if (iamp=='12'):
                axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC12']))
            if (iamp=='21'):
                axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC21']))  
            if (iamp=='22'):
                axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC22']))  
            x_start, x_end, y_start, y_end =\
                                             int(axis[2])-1, int(axis[3]), int(axis[0])-1, int(axis[1])
        if(Trimmed==True):
            x_start, x_end, y_start, y_end = self.OutRegion(iamp)
            logging.debug(str('trimmed = '+ str(Trimmed) +' '+ 'IlluRegion Amp ='+iamp+\
                    ' (x_start, x_end, y_start, y_end) = '+\
                    str(y_start)+' '+ str(y_end)+' '+ str(x_start)+' '+str(x_end)))
         
        return x_start, x_end, y_start, y_end


    def IsInAmp(self, x_coord, y_coord, Trimmed=False, **kwargs):
        for amp in self.amp:
            x_start, x_end, y_start, y_end = self.TelAmpLimits(amp, Trimmed=Trimmed)
            if((x_coord>=x_start) & (x_coord<x_end)\
               & (y_coord>=y_start) & (y_coord < y_end)):
                logging.debug(str('Coordinate (x,y) : '+str(y_coord)+' '+str(x_coord)+' '+' is in amp '+str(amp)))
                return amp
        return sys.exit('Coordinates outside of image')

                
    def Image(self, fitsimage):
        if(self.instrument=='ctio.9'):
            image_data = (pf.open(fitsimage))[0].data
            return image_data#[0:self.header['NAXIS2'],0 : self.header['NAXIS1']]

        
    def IlluRegion(self, frame, iamp) :
        x_start, x_end, y_start, y_end = self.TelAmpLimits(iamp=iamp)
        logging.debug(str('IlluRegion Amp ='+str(iamp)+' (x_start, x_end, y_start, y_end) = '+ str(y_start) +' '+ str(y_end)+' '+ str(x_start) +' '+ str(x_end)))
        img = frame[x_start:x_end, y_start:y_end]
        return frame[x_start:x_end, y_start:y_end]
       
          
    def OverscanRegion(self, frame, iamp) :
        if (iamp=='11'):
            name = 'BSEC11'
        if (iamp=='12'):
            name = 'BSEC12'
        if (iamp=='21'):
            name = 'BSEC21'
        if (iamp=='22'):
            name = 'BSEC22'
        axis = re.split(':|,', re.sub('\[|\]','', self.header[name]))
        logging.debug(str('OverscanRegion Amp ='+str(iamp)+' (x_start, x_end, y_start, y_end) = '+ str(int(axis[0])-1) +' '+ str(int(axis[1]))+' '+ str(int(axis[2])-1) +' '+ str(int(axis[3]))))  
        return frame[int(axis[2])-1 : int(axis[3]), int(axis[0])-1 : int(axis[1])]

    # Badly written, if other config occurs in future, add TSEC as self           
    def OutRegion(self, iamp):
        if (self.camera == 'single_amp'):
            axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC21']))
            sizeX21 = int(axis[3]) - (int(axis[2])-1)
            sizeY21 = int(axis[1]) - (int(axis[0])-1)
            out = [0, sizeX21, 0, sizeY21]
        else:
            axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC11']))
            sizeX11 = int(axis[3]) - (int(axis[2])-1)
            sizeY11 = int(axis[1]) - (int(axis[0])-1) 
            axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC12']))
            sizeX12 = int(axis[3]) - (int(axis[2])-1)
            sizeY12 = int(axis[1]) - (int(axis[0])-1) 
            axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC21']))
            sizeX21 = int(axis[3]) - (int(axis[2])-1)
            sizeY21 = int(axis[1]) - (int(axis[0])-1)
            axis = re.split(':|,', re.sub('\[|\]','',self.header['TSEC22']))
            sizeX22 = int(axis[3]) - (int(axis[2])-1)
            sizeY22 = int(axis[1]) - (int(axis[0])-1) 
            if (iamp=='11'):
                out = [0, sizeX11, 0, sizeY11]
            if (iamp=='12'):
                out = [0, sizeX12, sizeY11, sizeY11+sizeY12]
            if (iamp=='21'):
                out = [sizeX11, sizeX11+sizeX21, 0, sizeY21]
            if (iamp=='22'):
                out = [sizeX12, sizeX12+sizeX22, sizeY21, sizeY21+sizeY22]
        return out

    def XYoutSize(self):
        if (self.camera == 'single_amp'):
            self.XoutSize = 1024
            self.YoutSize = 1024
        else :
            self.XoutSize = 2048
            self.YoutSize = 2048
        return 
    
    def trim(self, frame):
        #outimg = np.zeros([self.lenY, self.lenX])
        outimg = np.zeros([self.YoutSize, self.XoutSize])
        for amp in self.amp:
            list = self.OutRegion(amp)
            outimg[list[0]:list[1], list[2]:list[3]] = self.IlluRegion(frame, amp)
        return outimg

    def OverscanSubtract_andTrim(self, frame):
        if ((str(self.header['COMMENT'])).find("Image is trimmed") >=0):
            sys.exit('Image is already trimmed, abort')
        #outimg = np.zeros([self.lenY, self.lenX])
        outimg = np.zeros([self.YoutSize, self.XoutSize])
        for amp in self.amp:
            list = self.OutRegion(amp)
            overscan = np.median(self.OverscanRegion(frame, amp))
            img      = self.IlluRegion(frame, amp)
            img  = img - overscan
            outimg[list[0]:list[1], list[2]:list[3]] = img
        return outimg


    '''From ds9 pix_x pix_y, return i, j in trimmed image, starting at 0,0'''
    '''Careful : x <-> y permutted ds9 <-> numpy.array'''
    def rawDS9toIJstart0(self, x_object, y_object):
        amp_size = 1024 # careful if new instruments !
        amp = self.IsInAmp( y_object-1, x_object-1, Trimmed=False)
        y0_amp, yf_amp, x0_amp, xf_amp\
            = self.TelAmpLimits(amp, Trimmed=False)

        out_x =np.floor( x_object /amp_size)*amp_size + x_object - x0_amp
        out_y =np.floor( y_object /amp_size)*amp_size + y_object - y0_amp
        return out_x, out_y
