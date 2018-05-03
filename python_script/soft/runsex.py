#!/usr/bin/env python 
'''
run sextractror on an image

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys
import toolbox as tb


def usage():
    print "run SExtractor on image and put it in outdir"
    print "Usage: runsex.py [outdir] [images]"
    print 
    sys.exit(1)

  

        
if __name__ == "__main__":
  
    if len(sys.argv) <= 1:
        usage()

    dir = sys.argv[1]
    image  = sys.argv[2:]
  
    for im in image:
        name   = im.split("/")[-1].replace(".fits","")
        outdir = os.path.join(dir,'prod',name)
        tb.RunSEx(im, outdir)
