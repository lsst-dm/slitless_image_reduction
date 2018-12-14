#!/usr/bin/env python
'''
remove flux at f(w) = f(w)-f(w/2)*.01
Author : Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import scipy.interpolate as interp


if __name__ == "__main__":
    show_stuff = True
    files = sys.argv[1:]
    wmin = 650.
    wmax = 1400.
    coef = 0.01
    for file in files:
        dict, values, names = tb.readcat(file)

        contaminant = np.array([values.field('w')*2, values.field('psf_GM_flux')*coef]).transpose()
        contaminant = contaminant[(contaminant[:, 0] >= wmin) & (contaminant[:, 0] <= wmax)]
        interp = interp.griddata(contaminant[:, 0], contaminant[:, 1], values.field('w'))
        interp[np.isnan(interp)] = 0
        corrected_flux = values.field('psf_GM_flux') - interp

        print('keys ', values.dtype.fields.keys())
        print('values ', values.item(1))
        print('values ', values.shape)
        outkey = names
        outkey.append('flux_decontaminated')

        out = []
        val = [list(i) for i in values]
        for i, j in zip(val, list(corrected_flux)):
            i.append(j)
            out.append(i)

        kk = list(dict.keys())
        kv = list(dict.values())

        outdict = {}
        for i, j in zip(kk, kv):
            outdict[i] = j[0]

        tb.DumpFile(outdict,
                    outkey,
                    zip(*out),
                    'test.list')
        print('writing : test.list')

        if show_stuff:
            print('values : ', dict.values())
            print('keys : ', dict.keys())
            #print values.field('w')
            print('Dict : ', dict.items())
            #print 'Array : ', values.field
            pl.plot(values.field('w'), values.field('psf_GM_flux'), label='raw flux')
            pl.plot(values.field('w'), corrected_flux, label='corrected')
            pl.legend()
            pl.show()
