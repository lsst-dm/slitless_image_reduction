#!/usr/bin/env python 
'''
return a trim, overscan subtracted image from a list of raw images

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os, sys
import numpy as np
import astropy.io.fits as pf
import telinst as instru
   

def usage():
    print "return a trim, overscan subtracted image from a list of raw images"
    print "Usage: median_stack.py [result_image_name] [list_of_images]"
    print 
    sys.exit(1)


 
if __name__ == "__main__":
    if len(sys.argv) <= 1:
        usage()
    OutIm = sys.argv[1]
    SOutIm = OutIm.replace('.fits', '_rms.fits')
    inputlist = sys.argv[2:]
    datas = []
    for img in inputlist :
        print 'opening ', img
        inst = instru.telinst(img)
        data = inst.Image(img)
        outimg = inst.OverscanSubtract_andTrim(data)
        datas.append(outimg)

    stack   = np.median(datas, axis=0)
    sstack  = np.std(datas, axis=0)
    hdr     = ((pf.open(inputlist[0]))[0].header).copy()
    info    = "Image is a stack of :"+ str(sys.argv[2:]) 
    hdr.add_comment(info)   
    pf.writeto(OutIm, stack, hdr, clobber=True)
    pf.writeto(SOutIm, sstack, hdr, clobber=True)
  
