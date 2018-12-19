#!/usr/bin/env python
'''
Test code for spectra extraction of slitless data

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''
from __future__ import print_function
import astropy.io.fits as pf
import numpy as np
import os
import math
import pylab as pl
import scipy
from scipy import optimize
import argparse
import copy
import sys
from croaks import NTuple
from matplotlib.colors import LogNorm

import lsst.utils
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import pysex


def result(fits):

    cat = pysex.run(fits, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_APER'],
                    conf_args={'PHOT_APERTURES': 5})
    print(cat['FLUX_APER'])

    image = afwImage.MaskedImageF(fits)

    binsize = 128
    nx = int(image.getWidth()/binsize) + 1
    ny = int(image.getHeight()/binsize) + 1
    bctrl = afwMath.BackgroundControl(nx, ny)

    bkgd = afwMath.makeBackground(image, bctrl)

    statsImage = afwMath.cast_BackgroundMI(bkgd).getStatsImage()

    image -= bkgd.getImageF(afwMath.Interpolate.NATURAL_SPLINE)

    return bkgd

    return


def grabargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pmax", type=int,
                        help="maximum p moment for the analysis",
                        default=3)
    args = parser.parse_args()
    return args


def main():
    #args = grabargs()
    inputfiles = sys.argv[1:]

    for image in inputfiles:
        result(image)


if __name__ == "__main__":

    main()
