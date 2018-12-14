#!/usr/bin/env python
'''
just this :
http://docs.astropy.org/en/stable/wcs/

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function
#from __future__ import division, print_function

from builtins import str
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
import telinst as instru

from astropy import wcs
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "pixel to sky"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract spectrum from raw image")
    parser.add_argument('-v', "--verbose",
                        help="verbose",
                        action='store_true')
    parser.add_argument('-i', "--reps", nargs='+', type=str,
                        help="list of reps to be processed",
                        default=None)
    args = parser.parse_args()
    return args


def runAstrometry(path, image):
    target = os.path.join(path, image)
    cmd = "solve-field %s --downsample 2 --scale-units arcminwidth --scale-low 10. --scale-high 15.0 --fits-image --overwrite --use-sextractor"%(
        target)
    print(cmd)
    os.system(cmd)
    remove = str(path+'calibrated-indx.png '+path+'calibrated-ngc.png '+path+'calibrated-objs.png '+path+'calibrated.axy '+path+'calibrated.corr ' +
                 path+'calibrated.match '+path+'calibrated.rdls '+path+'calibrated.solved '+path+'calibrated-indx.xyls ')+path+'calibrated.new '
    print('removing : ', remove)
    os.system('rm -f %s' % (remove))
    return


def copywcs(ref, target):
    hdulist1 = fits.open(ref)
    w = wcs.WCS(hdulist1[0].header)
    header = w.to_header()
    hdulist2 = fits.open(target, mode='update')
    target_header = hdulist2[0].header
    target_header.extend(list(header.items()))
    hdulist2.flush()
    return


if __name__ == "__main__":
    args = grabargs()
    reps = args.reps
    verbose = args.verbose

    if(verbose):
        Level = logging.getLogger().setLevel(logging.DEBUG)
        logging.debug('DEBUG mode')
    else:
        Level = logging.getLogger().setLevel(logging.INFO)
        logging.info('INFO mode')

    for rep in reps:
        image = 'calibrated.fits'
        target = os.path.join(rep, image)

        '''run astrometric match'''
        ''' write wcs to calibrated.fits'''
        runAstrometry(rep, image)
        ref = os.path.join(rep, 'calibrated.wcs')
        if os.path.isfile(ref):
            copywcs(ref, target)
        os.system('rm -f %s' % (ref))
