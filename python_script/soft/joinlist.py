#!/usr/bin/env python
'''
divide QE by CBP transmission 

Author: Augustin Guyonnet
aguyonnet@fas.harvardedu
'''
from __future__ import print_function

import os
import sys
import pylab as pl
import toolbox as tb
import numpy as np
from croaks import NTuple


if __name__ == "__main__":
    narg = len(sys.argv)
    if narg < 2:
        print("joinlist.py file1 file2")
        print("join points")
        print()

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    data1 = tb.readtxt(file1, ['img', 'wght', 'FLUX_AUTO'])
    data2 = tb.readtxt(file2, ['FILENAME', 'WAVELENGTH', 'SPOT_0_FLUX'])
    data1 = np.array(data1).transpose()
    data2 = np.array(data2).transpose()
    #data1 = data1[data1[:, 0].argsort()]
    #data2 = data2[data2[:, 0].argsort()]

    #parsing img number
    for i, item2 in enumerate(data2):
        img = ((item2[0].split("\\")[-1]).strip('.fits').split('_'))[-1]
        data2[i, 0] = img
        print(i, data2[i, 0])

    #sys.exit()
    img = []
    wght = []
    res = []

    for item1 in data1:
        for item2 in data2:
            if(int(item1[0]) == int(item2[0])):
                print('match ', item1[0], item1[1], item1[2], item2[1], item2[2])
                img.append(item1[0])
                wght.append(item1[1])
                res.append(float(item2[2])/float(item1[2]))

    pl.figure()
    zen = pl.figure()
    pl.plot(wght, res, 'r^')
    pl.ylim(3.02, 3.15)
    pl.xlabel('wavelength (nm)')
    pl.ylabel('Flux_Nick / Flux_Sextractor')
    pl.title('0th order light from CBP+Ronchi : Comparing two flux estimator')
    #pl.legend(loc='upper left')
    zen.savefig("cbp_ronchi_flux0th.pdf")
    pl.show()
