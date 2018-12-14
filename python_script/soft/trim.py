#!/usr/bin/env python
'''
trim illu section in fitsimage

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os
import sys
import astropy.io.fits as pf
import glob
import pylab as pl
import telinst as instru
#/Users/lpnhe/harvard/soft/python_script/build/lib/lsst/utils/


if __name__ == "__main__":

    Image = sys.argv[1:]
    for img in Image:
        inst = instru.telinst(img)
        data = inst.Image(img)
        outimg = inst.trim(data)
        hdr = (inst.header).copy()
        hdr.add_comment("Image is trimmed")
        (filepath, filename) = os.path.split(img)
        pf.writeto(os.path.join(inst.prod_path, 'trim_' + filename), outimg, hdr, clobber=True)
