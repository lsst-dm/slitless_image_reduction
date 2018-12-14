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


def extract_table(footprint, x_start, x_end, y_start, y_end, saturation):
    X = []
    Y = []
    end = x_end - x_start
    for i in range(end):

        remove_saturation = footprint[:, i]
        remove_saturation = remove_saturation[remove_saturation < saturation]

        flux = np.sum(remove_saturation)
        if flux > 0:
            X.append(i)
            Y.append(flux)
        else:
            print 'line ', i, ' is > 0.9* saturation'
    return X, Y


def graph(X, Y, Color='r', Linestyle='-', Name='raw.pdf', **kwargs):
    figOp = pl.figure(100)
    # Fait le graphique des largeurs en X VS position X pour la lecture a gauche
    g1 = pl.plot(X,
                 Y,
                 color=Color,
                 linestyle=Linestyle,
                 markersize=6,
                 label="Mes.")

    pl.xlabel('position (pixel)')
    pl.ylabel('flux (ADU)')
    pl.title('Raw spectra')

    figOp.savefig(Name)
    pl.clf()

    return 0


def result(image):
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
    x_start = 1400
    x_end = 2101
    y_start = 1140
    y_end = 1191
    footprint = image_data[y_start:y_end, x_start:x_end]

    # extract raw profile :
    X, Y = extract_table(footprint, x_start, x_end, y_start, y_end, saturation)

    # Draw raw profile :
    graph(X, Y, Color='r', Linestyle='-', Name=name + '.pdf')

    # Save values as list file :
    info = np.rec.fromrecords([([i1] + [i2])
                               for i1, i2,
                               in zip(X, Y)],
                              names=['X', 'Y'])
    info.view(NTuple).totxt(name + '.listj')

    #fig0=pl.figure(0)
    #pl.imshow(image_data, cmap='gray', norm=LogNorm())
    #cbar = pl.colorbar(ticks=[1.e3,1.e4,1.e5])
    #cbar.ax.set_yticklabels(['1,000','10,000','100,000'])
    #fig1=pl.figure(1)
    #pl.imshow(footprint, cmap='hot')
    ##pl.colorbar()

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
    pl.show()


if __name__ == "__main__":

    main()
