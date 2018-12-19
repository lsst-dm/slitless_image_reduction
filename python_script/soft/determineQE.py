#!/usr/bin/env python
'''
divide QE by CBP transmission 

Author: Augustin Guyonnet
aguyonnet@fas.harvardedu
'''

import os
import sys
import pylab as pl
import toolbox as tb
import numpy as np
import scipy.interpolate as interp

if __name__ == "__main__":

    step = 1.
    wght = np.arange(450, 901, step)

    qe = np.loadtxt(sys.argv[1])
    cbp = np.loadtxt(sys.argv[2])

    QE = interp.griddata(qe[:, 0], qe[:, 1], wght)
    CBP = interp.griddata(cbp[:, 0], cbp[:, 1], wght)

    out = QE/CBP

    filename = 'Oct2017QE.list'
    outdir = os.path.dirname(sys.argv[1])
    tb.DumpTuple(('wght', 'qe'), (wght, out), os.path.join(outdir, filename))

    pl.figure()
    pl.plot(qe[:, 0], qe[:, 1], 'r^', label='file1')
    pl.legend()

    pl.figure()
    pl.plot(cbp[:, 0], cbp[:, 1], 'bo', label='file2')
    pl.legend()

    pl.figure()
    pl.plot(wght, out, 'r^', label='cbp Norm')
    #pl.plot(wght, CBP, b'o')
    pl.legend()
    pl.show()
