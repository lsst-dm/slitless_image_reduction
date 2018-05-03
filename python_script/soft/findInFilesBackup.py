#!/usr/bin/env python 
'''
divide QE by CBP transmission 

Author: Augustin Guyonnet
aguyonnet@fas.harvardedu
'''

import os, sys
import pylab as pl
import toolbox as tb
import numpy as np
from croaks import NTuple


if __name__ == "__main__":
    narg = len(sys.argv)
    if narg<2 :
        print "findInFiles.py [file(s)] -v [value(s)]"
        print "return files if keywords are in it"
        print
    values = []
    Files   = []
    k = 1
    while( k<narg ):
        if( sys.argv[k][0] != "-" ):
            Files.append( sys.argv[k] )
            k += 1
        elif( sys.argv[k] == "-v" ):
            k += 1
            values=sys.argv[k:] 
            break
 
    nb = len(values)
    out = []
    for file in Files:
        found = 0
        data =  NTuple.fromtxt(file)
        for k in values:
            for key, val in data.keys.iteritems():                
                if ((str(val).lower()).find(k.lower())>=0):
                    found +=1
                    break
            if (found == 0):
                break
            
        if (found == nb):
            print file
            out.append(file)

    print 'files_list : ', out
   

    for spectrum in out:
        print spectrum
        data = NTuple.fromtxt(spectrum)
        fields = data.dtype.fields.keys()
        print fields
        x = np.array(data[:]['w'])
        y = np.array(data[:]['aper_flux'])
        pl.figure()
        pl.plot(x,y, 'r^')
        pl.show()
