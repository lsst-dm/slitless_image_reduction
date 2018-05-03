#!/usr/bin/env python 
'''
from a raw fitsimage return a spectrum

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys, re
import numpy as np
import pylab as pl
import logging
import toolbox as tb
import argparse
import logging
import reduceimage as ri




def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "process slitless images"
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract spectrum from raw image")
    parser.add_argument('-e',"--order", type=str, 
		        help = "either [m+1, m-1, all:default]",
                        default='all')
    parser.add_argument('-t',"--method", type=int, 
		        help = "either 1 or 2",
                        default=1)
    parser.add_argument('-r',"--raw", nargs='+', type=str, 
		        help = "raw fits image and x,y star position", 
		        default=None)
    parser.add_argument('-b',"--obj", type=str, 
		        help = "object name if not in fits header", 
		        default=None)
    parser.add_argument('-c',"--cal", 
		        help = "calibration using a synthetic flat", 
		        action='store_true')
    parser.add_argument('-k',"--dark", 
		        help = "subtract a masterdark", 
		        action='store_true')
    parser.add_argument('-d',"--disp",
		        help = "apply pixel to wavelength transformation", 
		      action='store_true')
    parser.add_argument('-m',"--map",
		        help = "build a map of defective pixels", 
		        action='store_true')
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
    parser.add_argument('-v',"--verbose", 
		        help = "verbose", 
		        action='store_true')
    parser.add_argument('-o',"--orient", type=str, 
		        help = "indicate dispersion direction", 
		        default='x')
    parser.add_argument('-i',"--imgs", nargs='+', type=str, 
	                help = "list of fitsimages to be processed", 
	                default=None)
    parser.add_argument('-s',"--isr",
		        help = "removal of periodic electronic pattern", 
		        action='store_true')
    parser.add_argument('-a',"--cosmic",
		        help = "replacing cosmics by median of surrounding pixels", 
		        action='store_true')
    parser.add_argument('-x',"--suffix", type=str, 
		        help = "appended to output rep name",
                        default='')
    parser.add_argument('-f',"--aperture", type=int, 
		        help = "aperture width",
                        default=150)
    args = parser.parse_args()
    return args


    

if __name__ == "__main__":
    args  = grabargs()
    if(args.raw):
        find_position = False
        x_object      = float(args.raw[0])
        y_object      = float(args.raw[1])
        Images        = [args.raw[2]]
    else :
        find_position =  True
        Images = args.imgs
    order         = args.order
    aperture      = args.aperture
    suffix        = args.suffix
    dark          = args.dark
    dispersion    = args.disp
    plot          = args.plot
    Map           = args.map
    object_name   = args.obj
    calibrate     = args.cal
    orientation   = args.orient
    verbose       = args.verbose
    method        = args.method
    isr           = args.isr
    cosmic        = args.cosmic
    if(verbose):
        Level=logging.getLogger().setLevel(logging.DEBUG)
        logging.debug('DEBUG mode')
    else:
        Level=logging.getLogger().setLevel(logging.INFO)
        logging.info('INFO mode')
     

    for Image in Images:
        ''' Overscan subtract, trim and write new image in prod rep '''
        ''' Subtract master-dark -- scaled to same exposure time    '''
        logging.info('Processing : '+ Image)
        prepare = ri.Prepare(Image,
                             object_name = object_name,
                             suffix = suffix)
        if prepare.Check():
            continue
        
        prepare.TrimOverscanSubtract()
        prepare.Dark(dark)
        prepare.RunSEx(os.path.join(prepare.outdir, 'calibrated.fits'))
        logging.debug('OUTDIR ' + prepare.outdir)
        '''if object pos not provided, find it in se.list 
        from Sextractor on calibrated.fits  '''
        if (find_position == True):
            x_object, y_object = prepare.FindObject()
        else:
            '''From ds9 pixel i pixel j, return i, j in trimmed image, starting at 0,0'''
            x_object, y_object = prepare.instrument.rawDS9toIJstart0(x_object, y_object)
        ''' run sectractor on image to extract a background map  '''
        ''' handling the extraction of the spectrum              '''

        ### On CTIO slitless images, Seeing does not find stars but hot pixels
        ### --> hacked on default value 1.5
        #seeing, sx, sy  = tb.Seeing(os.path.join(prepare.outdir, 'se.list'), plot)
        seeing = 1.5 ; sx = 0. ; sy = 0. 
        process = ri.Process(Image,
                             Map         = Map,
                             orientation = orientation,
                             dispersion  = dispersion,
                             plot        = plot,
                             method      = method,
                             calibrate   = calibrate,
                             outdir      = prepare.outdir,
                             cosmic      = cosmic,
                             isr         = isr,
                             aperture    = aperture,
                             object_name = object_name)
              
     

        if (order=='all'):
            orders = ['m+1', 'm-1']
        else:
            orders = [order]
  
        pixel, wavelength, raw_profile, calib_profile  =\
                  process.SpectrumExtraction(x_object, y_object,
                                             seeing=seeing,
                                             sx=sx, sy=sy,
                                             orders = orders,
                                             offset = 100)
 
        
        if plot:
            print 'show plots'
            pl.show()
