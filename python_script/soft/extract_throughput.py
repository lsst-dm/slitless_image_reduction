#!/usr/bin/env python 
'''
from monochromatic fitsimages and monitoring photodiode files, 
return a QE*optical throughput curve

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys
import astropy.io.fits as pf
import glob
import pylab as pl
import numpy as np
from scipy import stats
import telinst as instru
import croaks
from croaks import NTuple



def usage():
    print 'extract_throughput.py [file]'
    print "file is a list of fitsimages and their associated Iph files"
    print "from monochromatic fitsimages and monitoring photodiode files," 
    print "return a QE*optical throughput curve"
    print "add [Plot] to see check plots: "
    print


def readfile(inputfile, dic=[]):
    data = croaks.NTuple.fromtxt(inputfile)
    value =[]
    for name in dic:
        value.append(np.array(data[:][name]))
    return value


def initialize(li):
    ronchi  = 0
    fits    = li[0] # image
    file    = li[1] # img photodiode monitoring
    fileD   = li[3] # dark photodiode monitoring
    inst    = instru.telinst(fits)
    data    = inst.Image(fits)
    img     = inst.header['IMAGETYP']
    object  = inst.header['OBJECT']
    exptime = inst.header['EXPTIME']
    filters = inst.header['FILTERS']
    return ronchi, fits, file, inst, data, img, object, exptime, filters, fileD


def SubtractDark(data, exptime, darkfilename):
    darkframe   = instru.telinst(darkfilename) # the dark image
    exptimedark = darkframe.header['EXPTIME']
    scaling     = exptimedark / exptime 
    arr         = darkframe.Image(darkfilename)
    data       -= (arr * scaling)
    print 'dark :', darkfilename, ' scaling factor = ', scaling
    return data


def Photodiode(iph, exptime):
    n    = len(iph[:,0])
    norm = float(exptime) / n
    Iph  = np.sum(iph[:,0]) * norm
    sIph = np.std(iph[:,0])/ np.sqrt(n) * float(exptime)
    return -Iph, sIph


def FluxInAmp(amp, data, plot=False, **kwargs):
    list = inst.OutRegion(amp)
    frame = data[list[0]:list[1], list[2]:list[3]]
    pix = (frame).flatten()
    fact = 5
    pix2, low, upp = stats.sigmaclip(pix, fact, fact)
    print 'len pix = ', len(pix), ', len pix2 = ', len(pix2)
    if(plot==True):
        pl.hist(pix, histtype='bar', log= True, range=(-4000, 4000), bins = 100, color='green',label = 'trial')
        pl.hist(pix2, histtype='bar', log= True, range=(-4000, 4000), bins = 100, color='crimson',label = 'trial')
        pl.show()

    mean  = np.median(pix2) 
    smean = np.std(pix2)

    '''Flag clipped pixels'''

    np.clip(frame, low, upp, out = data[list[0]:list[1], list[2]:list[3]])
    print 'minmax outimage : ', np.amax(data[list[0]:list[1], list[2]:list[3]]), np.amin(data[list[0]:list[1], list[2]:list[3]])
    
    return mean, smean, low, upp


def DumpFile(string, name):
    info = np.rec.fromrecords([([i1] +[i2]+[i3]+[i4]+[i5]+[i6]+[i7]+[i8])
                               for i1, i2, i3, i4, i5, i6, i7, i8
                               in string],
                              names = ['wavelength', 'ronchi', 'amp', 'mean', 'smean', 'IphTot', 'sIphTot', 'exptime'])
    info.view(NTuple).totxt(name+ '.list' )


'''Read Img photodiode file'''
def Iph(file):   
    with open(file, 'r') as f:
        comments = f.readline()
    iph        = np.loadtxt(file)
    comments   = comments.split() 
    wavelength = comments[comments.index('WAVELENGTH:')+1]
    '''Integrated photocurrent'''
    ph, sph = Photodiode(iph, exptime)
    return ph, sph, wavelength


    

if __name__ == "__main__":
    dark      = False
    name      = 'calibADU'
    result    = []    
    inputfile = sys.argv[1]
    Plot      = False
    if (str(sys.argv[1:]).find('Plot') >= 0):
        Plot = True
    WriteFits = False
    if (str(sys.argv[1:]).find('WriteFits') >= 0):
        WriteFits = True
    
 
    if (str(sys.argv[1:]).find('Dark') >= 0):
        dark = True
        list = readfile(inputfile, ['image','photodiode', 'dark', 'dark_pd'])
    else:
        list = readfile(inputfile, ['image','photodiode'])


        
    for li in zip(*list):
        ronchi, fits, file, inst, data, img, object, exptime, filters, fileD = initialize(li)
        print 'Processing : ' , fits
        print filters
        if not (filters.find('ronchi')): # to distinguish datasets
            print 'ronchi in path -> 1'
            ronchi =1  # if ronchi is in 
    
        '''Subtract dark'''
        if (dark==True):
            print 'dark=true'
            data = SubtractDark(data, exptime, li[2])

          
        ''' get Amp*exptime from Img and its associated dark'''
        iph, siph, wavelength    = Iph(file)
        Diph, Dsiph, Dwavelength = Iph(fileD)
        
        print 'Img and Dark Photocurrents : ',iph, siph, wavelength, Diph, Dsiph, Dwavelength

        '''Subtract dark current '''
        iph = iph - Diph
        siph= np.sqrt(siph*siph + Dsiph*Dsiph)
        print 'Iph dark subtracted : ', iph, siph
        
        if(WriteFits==True):
            hdr    = (inst.header).copy()  
            hdr.add_comment("Dark subtracted")
            hdr.set('LAMBDA', wavelength, 'Wavelength of illumination (nm)')
            hdr.set('IPHTOT',iph, 'Integrated Img ph - Integ.dark ph (A*s)')
            hdr.set('SIPHTOT',siph, 'rms/sqrt(n) of IPHTOT')

            
            '''Flux in amps'''
        for amp in inst.amp :
            mean, smean, low, high = FluxInAmp(amp, data, plot=Plot)
            result.append([wavelength, ronchi, amp, mean, smean, iph, siph, exptime] )
            print 'wavelength, ronchi, amp, mean, smean, Iph, sIph, exptime ', wavelength, ronchi, amp, mean, smean, iph, siph, exptime
            if(WriteFits==True):
                hdr.set('LOW'+amp,low, 'Pix min value after 5s clipping for the amp')
                hdr.set('HIGH'+amp,high, 'Pix max value after 5s clipping for the amp')

            
 
        '''write fitsimage with clip low and high'''
        '''When using image for throughput, I should do something to these min/max pix'''
        if(WriteFits==True):
            pf.writeto('dark_' + fits
                       , data, hdr, overwrite=True)


    # Return calibration as file :
    DumpFile(result, name)
     
