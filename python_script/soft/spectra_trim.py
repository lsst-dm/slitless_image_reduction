#!/usr/bin/env python 
'''
Test code for spectra extraction of slitless data

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''
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



def extract_profile(footprint, mask, direction):
    Y =[]
    if(direction=="y"):
        footprint = footprint.T
        mask      = mask.T
    rows = len(footprint)
    cols = len(footprint[0])
    print 'rows = ', rows, ' cols = ', cols
    for i in range(0,rows): 
        keep = footprint[i , :]
        remo = mask[i , :]
        keep = keep[remo==0]

        flux = np.median(keep)  
        Y.append(flux)
    return Y
    

def extract_table(footprint, x_start, x_end, y_start, y_end, saturation, profile):
    X =[] ; Y =[]
    end = x_end - x_start
 
    for i in range(end):
        remove_saturation = footprint[: , i]
        
        if (len(profile)!= 0) :
            remove_saturation = remove_saturation - profile[i]
          
        remove_saturation = remove_saturation[remove_saturation<saturation]
        flux = np.sum(remove_saturation)
        if flux > 0 :
            X.append(i)      
            Y.append(flux)
        else :
            print 'line ', i , ' is > 0.9* saturation'
    return X, Y
    

    

def graph(X, Y, Color = 'r', Linestyle = '-' , Name = 'raw.pdf', **kwargs): 
    figOp=pl.figure(100)
    # Fait le graphique des largeurs en X VS position X pour la lecture a gauche 
    g1 = pl.plot(X, 
                 Y,
                 color= Color,
                 linestyle = Linestyle,
                 markersize = 6, 
                 label = "Mes.")

    pl.xlabel('position (pixel)')
    pl.ylabel('flux (ADU)')
    pl.title('Raw spectra')
    
    figOp.savefig(Name)
    pl.show()
    pl.clf()
    
    return 0
    

def result(image, mask, direction):
    print ' Opening of image ', image
    f = pf.open(image)
    head = f[0].header
    #f.info()
    name = head['OBJECT'] + head['RECID'] 
    image_data = f[0].data
    f.close()
    

    # get saturation :
    saturation = np.max(image_data)*0.9
    # Select the footprint :
    x_start = 1250
    x_end   = 2001
    y_start = 1130
    y_end   = 1200 
    footprint = image_data[y_start:y_end, x_start:x_end]
    profile = []
    if(mask):
        name += '_pro'
        m = (pf.open(mask))[0].data
        m = m[1030:2000, x_start:x_end]
        f = image_data[1030:2000, x_start:x_end]
        profile = extract_profile(f, m, direction)
        
    # extract raw profile :
    print len(profile)
    X, Y = extract_table(footprint, x_start, x_end, y_start, y_end, saturation, profile)

    # Draw raw profile :   
    graph(X, Y, Color = 'r', Linestyle = '-' , Name = name+ '.pdf' )

    # Save values as list file :
    info = np.rec.fromrecords([([i1] + [i2])
                           for i1, i2, 
                               in zip( X, Y)],
                          names = ['X', 'Y'])
    info.view(NTuple).totxt(name+ '.list' )

    
    #fig0=pl.figure(0)
    #pl.imshow(image_data, cmap='gray', norm=LogNorm())
    #cbar = pl.colorbar(ticks=[1.e3,1.e4,1.e5])
    #cbar.ax.set_yticklabels(['1,000','10,000','100,000'])
    #fig1=pl.figure(1)
    #pl.imshow(footprint, cmap='hot')
    ##pl.colorbar()
    
    return 



    
def main():
    inputfiles = []
    mask = None
    direction = None
    narg = len(sys.argv)
    if narg<2 :
        print "header.py [fitsimage(s)] -k [mask] [direction : 'x' or 'y']"
        print "If k is none, no mask is used"
        print
    keywords = []
    Images   = []
    k = 1
    while( k<narg ):
        if( sys.argv[k][0] != "-" ):
            inputfiles.append( sys.argv[k] )
            k += 1
        elif( sys.argv[k] == "-k" ):
            k += 1
            mask = sys.argv[k]
            direction = sys.argv[k+1]
            break
    
    for image in inputfiles :
        result(image, mask, direction)
        


    
if __name__ == "__main__":

	main()
