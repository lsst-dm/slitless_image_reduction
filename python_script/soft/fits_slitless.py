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

from matplotlib.colors import LogNorm

def extract_table(footprint, x_start, x_end, y_start, y_end, saturation):
    X =[] ; Y =[]
    end = x_end - x_start
    for i in range(end):
   
        remove_saturation = footprint[: , i]
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
    
    return 0
    

def result(image, xini , xfin, yini, yfin):
    print ' Opening of image ', image
    f = pf.open(image)
    head = f[0].header
    detsec22 = head['ADSEC22']
    print 'detsec22', detsec22

    #f.info()
    name = head['OBJECT'] + head['RECID'] + '.pdf'
    
    image_data = f[0].data
    print(image_data.shape)
    
   
    saturation = np.max(image_data)*0.9

    # Select the footprint :
    x_start = int(xini) # 1400
    x_end   = int(xfin) # 2101
    y_start = int(yini) # 1140
    y_end   = int(yfin) # 1191  
    footprint = image_data[y_start:y_end, x_start:x_end]

    # extract raw profile :
    X, Y = extract_table(footprint, x_start, x_end, y_start, y_end, saturation)
    # Draw raw profile :   
    graph(X, Y, Color = 'r', Linestyle = '-' , Name = name )

    
    fig0=pl.figure(0)
    pl.imshow(image_data, cmap='gray', norm=LogNorm())
    cbar = pl.colorbar(ticks=[1.e3,1.e4,1.e5])
    cbar.ax.set_yticklabels(['1,000','10,000','100,000'])

    fig1=pl.figure(1)
    #image_data = image_data[image_data>11]
    #pl.imshow(image_data, cmap='hot', clim = (10,40))
    pl.imshow(image_data, cmap='hot', vmin=5, vmax=10)
    pl.colorbar()
    pl.show()
    sys.exit()
    
    pix = (f[0].data).flatten()
    fig2=pl.figure(2)
    saturation = np.max(pix)*0.9
    pix = pix[pix<saturation]
    print 'number of pixels = ', len(pix)
    mea = 'mean = '+ str(np.mean(pix))
    print mea
    pl.hist(pix, log=True, range=(0, 40), histtype='bar', bins = 100, color='crimson',label = 'trial'); pl.show()
    pl.title('Master-dark residuals')
    pl.xlabel('count (ADU)')
    
    fig2.savefig("hist.pdf") 
    
    #pl.hist(pix, log = True, bins = 100); pl.show() 
    return 


def grabargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("--pmax", type=int, \
			    help = "maximum p moment for the analysis", \
			    default=3)
	args = parser.parse_args()
	return args



    
def main():
    #args = grabargs()
    xini = sys.argv[1]
    xfin = sys.argv[2]
    yini = sys.argv[3]
    yfin = sys.argv[4]
    inputfiles = sys.argv[5:]

    for image in inputfiles :
        result(image, xini, xfin, yini, yfin)

        


    
if __name__ == "__main__":

	main()
