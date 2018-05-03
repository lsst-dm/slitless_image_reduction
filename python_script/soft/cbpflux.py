#!/usr/bin/env python 
'''
run sextractror on an image

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys
import toolbox as tb
from croaks import NTuple
import astropy.io.fits as pf

def zeroth(cat):
    data = tb.ParseSexCat(cat)
    data = data[(data.field('X_IMAGE')>338) & (data.field('X_IMAGE')<352)]
    data = data[(data.field('Y_IMAGE')>590) & (data.field('Y_IMAGE')<610)]
    if(len(data)!=1):
        print 'WARNING : ', len(data), ' object found '
    return data


def usage():
    print "look at flux of spot in se.list files"
    print "Usage: cbpflux.py [files]"
    print 
    sys.exit(1)

  

        
if __name__ == "__main__":
  
    if len(sys.argv) <= 1:
        usage()

    files = sys.argv[1:]
    data = []
    names = zeroth(files[0]).dtype.names
    for file in files:
        dirname  = os.path.dirname(file)
        img      = dirname.split("_")[-1]
        fits     = os.path.join(dirname,'segmentation.fits')
        filters  = pf.open(fits)[0].header.get('FILTERS')
        comments = str(pf.open(fits)[0].header.get('COMMENT'))
        wght     = comments.split("WAVE")[1].split(" ")[1]
        print fits, 'wght = ', wght
        obj = zeroth(file)
        if (len(obj)!=1):
            continue
        obj = list(obj[0].tolist())
        obj = [int(img)] + [float(wght)] + obj
        data.append(obj)


    names = ['img']+['wght']+list(names)
    tb.DumpTuple(names, zip(*data), 'cbp.list')
