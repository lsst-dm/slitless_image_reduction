#!/usr/bin/env python
'''
return the median horizontal profile of a frame

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''
from __future__ import print_function

import os
import sys
import astropy.io.fits as pf
import numpy as np
import glob
import pylab as pl
import telinst as instru
#/Users/lpnhe/harvard/soft/python_script/build/lib/lsst/utils/


def graph(X, Y, Color='r', Linestyle='-', Name='raw', **kwargs):

    # Fait le graphique des largeurs en X VS position X pour la lecture a gauche
    g1 = pl.plot(X,
                 Y,
                 color=Color,
                 linestyle=Linestyle,
                 markersize=6,
                 label=Name)

    return 0


def extract_profile(footprint, mask, direction):
    X = []
    Y = []
    if(direction == "y"):
        footprint = footprint.T
        mask = mask.T
    rows = len(footprint)
    cols = len(footprint[0])
    print('rows = ', rows, ' cols = ', cols)
    for i in range(0, rows):
        keep = footprint[i, :]
        remo = mask[i, :]
        keep = keep[remo == 0]

        flux = np.median(keep)
        X.append(i)
        Y.append(flux)
    return X, Y


if __name__ == "__main__":

    direction = sys.argv[1]# 'x' is default, 'y' transpose the array

    image = sys.argv[2]
    mask = sys.argv[3]

    fig = pl.figure(1)
    inst = instru.telinst(image)
    data = inst.Image(image)
    exclude = inst.Image(mask)
    footprint = data[1030:2000, 1030:2000]
    mask = exclude[1030:2000, 1030:2000]
    X, Y = extract_profile(footprint, mask, direction)

    name = "bkgd_" + \
           (str(direction) + inst.header['RECID']).replace('.', '_')#+inst.header['IMAGETYP']

    # Draw raw profile :
    graph(X, Y, Color='r', Linestyle='-', Name=name)

    pl.xlabel('position (pixel)')
    pl.ylabel('counts (ADU)')
    #pl.title('Overscans')
    name += '.pdf'
    pl.show()
    fig.savefig(name)
