#!/usr/bin/env python
'''
Looking at a PanStarrs Image

'''
from __future__ import print_function

from builtins import str
from builtins import range
from builtins import object
import os
import sys
import astropy.io.fits as pf
import pylab as pl
import re
import numpy as np


class panstarrs(object):
    def __init__(self, fitsimage, Level=None, **kwargs):
        list = self.Acceptor(fitsimage)
        self.instrument = list[0]
        self.nccd = list[1]
        self.namp = list[2]
        self.lenX = list[3]
        self.lenY = list[4]
        self.fits = list[5]
        self.prod_path = list[6]
        self.amp = list[7]
        self.flip = False

    def Acceptor(self, fitsimage):
        fits = pf.open(fitsimage)
        head = fits[0].header
        if (str(head.get('TELESCOP') == 'PS1')):
            nccd = 1
            namp = int(head.get('NAMPS'))
            nextend = int(head.get('NEXTEND'))
            if(namp != nextend):
                sys.exit('namp != nextend')
            lenX = int(head['IMNAXIS1'])
            lenY = int(head['IMNAXIS2'])
            instrument = 'ps1'
            amp = list(range(1, namp+1))
            prod_path = './'
            return [instrument, nccd, namp, lenX, lenY, fits, prod_path, amp]
        else:
            return NULL

    def IlluRegion(self, frame, iamp):
        x_start, x_end, y_start, y_end = self.TelAmpLimits(iamp=iamp)
        return frame[x_start:x_end, y_start:y_end]

    def TelAmpLimits(self, iamp=None, Trimmed=False, **kwargs):
        if (Trimmed == False):
            axis = re.split(':|,', re.sub('\[|\]', '',
                                          self.fits[iamp].header.get('datasec')))
        else:
            sys.exit('Amp limits not defined for trimmed img yet')
        x_start, x_end, y_start, y_end =\
            int(axis[2])-1, int(axis[3]), int(axis[0])-1, int(axis[1])
        return x_start, x_end, y_start, y_end

    def OverscanRegion(self, frame, iamp):
        axis = re.split(':|,', re.sub('\[|\]', '',
                                      self.fits[iamp].header.get('biassec')))
        return frame[int(axis[2])-1: int(axis[3]), int(axis[0])-1: int(axis[1])]

    def OutRegion(self, iamp):
        self.flip = False
        axis = re.split(':|,', re.sub('\[|\]', '',
                                      self.fits[iamp].header.get('detsec')))
        x_start, x_end = int(axis[2])-1, int(axis[3])

        if(int(axis[0]) > int(axis[1])):

            self.flip = True
            y_start, y_end = int(axis[1])-1, int(axis[0])
        else:
            y_start, y_end = int(axis[0])-1, int(axis[1])
        return [x_start, x_end, y_start, y_end]

    def trim(self, gain=False, **kwargs):
        outimg = np.zeros([self.lenY, self.lenX])
        for amp in self.amp:
            list = self.OutRegion(amp)
            frame = self.IlluRegion(self.fits[amp].data, amp)
            print('Amp, Boundaries : ', amp, list)
            if self.flip == True:
                frame = np.fliplr(frame)
            if gain is True:
                print(amp, self.Gain(amp))
                frame *= self.Gain(amp)
            outimg[list[0]:list[1], list[2]:list[3]] = frame
        return outimg

    def Mask(self, **kwargs):
        outimg = np.zeros([self.lenY, self.lenX])
        for amp in self.amp:
            list = self.OutRegion(amp)
            x_start, x_end, y_start, y_end = self.TelAmpLimits(amp)
            frame = np.ones([x_end-x_start, y_end - y_start])
            outimg[list[0]:list[1], list[2]:list[3]] = frame
        return outimg

    def OverscanSubtract_andTrim(self, gain=False, **kwargs):
        if ((str(self.fits[0].header.get('COMMENT'))).find("Image is trimmed") >= 0):
            sys.exit('Image is already trimmed, abort')
        outimg = np.zeros([self.lenY, self.lenX])
        for amp in self.amp:
            list = self.OutRegion(amp)
            overscan = np.median(self.OverscanRegion(self.fits[amp].data, amp))
            img = self.IlluRegion(self.fits[amp].data, amp)
            img = img - overscan
            if gain is True:
                img *= self.Gain(amp)
            if self.flip == True:
                img = np.fliplr(img)
            outimg[list[0]:list[1], list[2]:list[3]] = img
        return outimg

    def Gain(self, amp):
        gain = self.fits[amp].header.get('GAIN')
        if gain is None:
            gain = 1.
            print('amp = ', amp, ' No gain found, default=1')
        return gain


if __name__ == "__main__":

    fits_file = 'c7675g0008o44.fits'
    ps1 = panstarrs(fits_file)
    #outimg    = ps1.trim(gain=True)
    outimg = ps1.OverscanSubtract_andTrim(gain=True)
    mask = ps1.Mask()
    hdr = (ps1.fits[0].header).copy()

    hdr.add_comment("Image is trimmed")
    hdr['IMG_NAME'] = (fits_file, 'initial image name')
    outname = str('trim_' + fits_file)
    pf.writeto(outname, outimg, hdr, clobber=True)

    outmask = str('mask_' + fits_file)
    comp = pf.CompImageHDU(name=outmask, data=mask, header=hdr, compression_type='PLIO_1')
#    comp = pf.CompImageHDU(name=outmask, data=mask, header=hdr, compression_type='RICE_1')
    comp.writeto(outmask, clobber=True)
