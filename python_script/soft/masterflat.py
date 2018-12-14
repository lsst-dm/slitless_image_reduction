#!/usr/bin/env python
'''
Take the median of a stack of normalized images

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os
import sys
import numpy as np
import astropy.io.fits as pf
import telinst as instru


def usage():
    print "return a masterflat from median combining normalized flats"
    print "Usage: masterflat.py [result_image_name] [masterbias_subtracted2eachflat] [list_of_flat]"
    print
    sys.exit(1)


if __name__ == "__main__":

    if len(sys.argv) <= 1:
        usage()

    outim = sys.argv[1]
    inputlist = sys.argv[3:]
    datas = []
    darkstack = (pf.open(sys.argv[2]))[0].data

    for img in inputlist:
        inst = instru.telinst(img)
        data = inst.Image(img)
        outimg = inst.OverscanSubtract_andTrim(data) - darkstack
        outimg = outimg / np.median(outimg)
        datas.append(outimg)

    stack = np.median(datas, axis=0)
    hdr = ((pf.open(inputlist[0]))[0].header).copy()
    info = "Image is a stack of :" + str(sys.argv[2:])
    hdr.add_comment(info)
    pf.writeto(outim, stack, hdr, overwrite=True)
