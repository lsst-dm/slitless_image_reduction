#!/usr/bin/env python 
'''
just this :
http://docs.astropy.org/en/stable/wcs/

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
#from __future__ import division, print_function

import os, sys, re
import numpy as np
import pylab as pl
import logging
import toolbox as tb
import argparse
import logging
import reduceimage as ri
import telinst as instru

from astropy import wcs
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import ICRS
from astropy.coordinates import match_coordinates_sky



def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "cpy and paste wcs from img -i to img -t"
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract spectrum from raw image")
    parser.add_argument('-i',"--ref", type=str, 
	                help = "reference wcs", 
	                default=None)
    parser.add_argument('-t',"--targets", nargs='+', type=str, 
		        help = "where to past it",
                        default='')
    args = parser.parse_args()
    return args


    

if __name__ == "__main__":
    args    = grabargs()
    targets = args.targets
    ref     = args.ref

    hdulist1 = fits.open(ref)
    w        = wcs.WCS(hdulist1[0].header)
    header   = w.to_header()
    print header.values(), header.keys(), header.items()

    for target in targets :
        hdulist2      = fits.open(target, mode='update')
        target_header = hdulist2[0].header
        target_header.extend(header.items())
        print(repr(header))
        hdulist2.flush()


