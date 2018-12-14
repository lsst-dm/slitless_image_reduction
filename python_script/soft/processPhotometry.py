#!/usr/bin/env python
'''
from a raw fitsimage return a spectrum

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os
import sys
import re
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

    parser.add_argument('-r', "--raw", nargs='+', type=str,
                        help="raw fits image and x,y star position",
                        default=None)
    parser.add_argument('-b', "--obj", type=str,
                        help="object name if not in fits header",
                        default=None)
    parser.add_argument('-c', "--cal",
                        help="calibration using a synthetic flat",
                        action='store_true')
    parser.add_argument('-k', "--dark",
                        help="subtract a masterdark",
                        action='store_true')
    parser.add_argument('-d', "--disp",
                        help="apply pixel to wavelength transformation",
                        action='store_true')
    parser.add_argument('-m', "--map",
                        help="build a map of defective pixels",
                        action='store_true')
    parser.add_argument('-p', "--plot",
                        help="show control plots",
                        action='store_true')
    parser.add_argument('-v', "--verbose",
                        help="verbose",
                        action='store_true')
    parser.add_argument('-i', "--imgs", nargs='+', type=str,
                        help="list of fitsimages to be processed",
                        default=None)
    parser.add_argument('-s', "--isr",
                        help="removal of periodic electronic pattern",
                        action='store_true')
    parser.add_argument('-a', "--cosmic",
                        help="replacing cosmics by median of surrounding pixels",
                        action='store_true')
    parser.add_argument('-x', "--suffix", type=str,
                        help="appended to output rep name",
                        default='')
    parser.add_argument('-g', "--proddir", type=str,
                        help="prod dir file",
                        default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = grabargs()
    if(args.raw):
        find_position = False
        x_object = float(args.raw[0])
        y_object = float(args.raw[1])
        Images = [args.raw[2]]
    else:
        find_position = True
        Images = args.imgs
    proddir = args.proddir
    suffix = args.suffix
    dark = args.dark
    dispersion = args.disp
    plot = args.plot
    Map = args.map
    object_name = args.obj
    calibrate = args.cal
    verbose = args.verbose
    isr = args.isr
    cosmic = args.cosmic
    if(verbose):
        Level = logging.getLogger().setLevel(logging.DEBUG)
        logging.debug('DEBUG mode')
    else:
        Level = logging.getLogger().setLevel(logging.INFO)
        logging.info('INFO mode')

    for Image in Images:
        ''' Overscan subtract, trim and write new image in prod rep '''
        ''' Subtract master-dark -- scaled to same exposure time    '''
        logging.info('Processing : ' + Image)
        prepare = ri.Prepare(Image,
                             object_name=object_name,
                             suffix=suffix,
                             proddir=proddir)
        if prepare.Check():
            continue

        prepare.TrimOverscanSubtract()
        prepare.Dark(dark)

        '''flatfielding'''
        name = prepare.Flat()
        print name
        #prepare.RunSEx(os.path.join(prepare.outdir, name))
        prepare.RunSEx(os.path.join(prepare.outdir, 'calibrated.fits'))
        logging.info('written in : ' + prepare.outdir)

        if plot:
            print 'show plots'
            pl.show()
