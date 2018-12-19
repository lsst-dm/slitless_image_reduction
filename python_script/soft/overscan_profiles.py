#!/usr/bin/env python
'''
return overscan profile

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

from builtins import str
from builtins import next
from builtins import range
import os
import sys
import astropy.io.fits as pf
import numpy as np
import glob
import pylab as pl
import telinst as instru
import matplotlib.cm as cm
import toolbox as tb
#/Users/lpnhe/harvard/soft/python_script/build/lib/lsst/utils/


def graph(X, Y, Color='r', Linestyle='-', Name='raw', **kwargs):

    # Fait le graphique des largeurs en X VS position X pour la lecture a gauche
    g1 = pl.plot(X,
                 Y,
                 color=Color,
                 linestyle=Linestyle,
                 markersize=4,
                 label=Name)

    return 0


def extract_profile(footprint, direction):
    X = []
    Y = []
    if(direction == "y"):
        footprint = footprint.T
    rows = len(footprint)
    cols = len(footprint[0])
    print('rows = ', rows, ' cols = ', cols)
    for i in range(0, rows):
        rem = footprint[i, :]
        flux = np.mean(rem)
        X.append(i)
        Y.append(flux)
    return X, Y


if __name__ == "__main__":
    amps = ('11', '12', '21', '22')
    direction = sys.argv[1]# 'x' is default, 'y' transpose the array
    Image = sys.argv[2:]
    for img in Image:
        figOp = pl.figure(1)
        inst = instru.telinst(img)
        data = inst.Image(img)
        (filepath, filename) = os.path.split(img)
        name = filename.split('.')[0]+str(direction)

        colors = iter(cm.rainbow(np.linspace(0, 1, len(amps))))
        for amp in amps:
            overscan = inst.OverscanRegion(data, amp)
            X, Y = extract_profile(overscan, direction)

            # Draw raw profile and dump tuple:
            graph(X, Y, Color=next(colors), Linestyle='-', Name=str(amp))
            tb.DumpTuple(['X', 'Y'], [X, Y], name+str(amp)+'.list')

        pl.xlabel('position (pixel)')
        pl.ylabel('counts (ADU)')
        pl.title('Overscans' + name)
        pl.legend(loc='upper left')
        name += '.pdf'
        figOp.savefig(name)
