#!/usr/bin/env python
'''
return plots of equivalent width

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function

from builtins import next
from builtins import str
import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import croaks
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import pandas as pd


def func(x, a, alpha):
    return a * x**alpha


if __name__ == "__main__":

    airmass = []
    jd = []
    wgth = []
    ew = []
    sew = []
    img = []
    seeing = []
    zenith = []
    files = sys.argv[1:]

    reference = []
    colors = iter(cm.rainbow(np.linspace(0, 1, len(files)+1)))
    fig, ax1 = pl.subplots(figsize=(7, 5))
    #pl.xlim([350,1000])
    pl.xlim([700, 1000])
    pl.ylim([0.9, 1.1])
    pl.ylabel('aper_flux / aper_100')
    pl.xlabel('wght (nm)')

    for file in files:
        values, keys = tb.readlist(file, ['pixel', 'w', 'aper_flux'])
        filters = keys['FILTERS']
        obj = re.split('[ /]', keys['IMG_NAME'])
        img = obj[-2]
        print()
        print('raw file : ', obj[0], img, filters)
        print(file)
        aperture = file.split('/')[-2].split('_')[-1].split('aper')[-1]
        print(aperture)
        if (str(aperture) == '100'):
            print('reference is aper 100')
            reference = values[2]
            wght = values[0]

        print(len(wght), len(values[0]))
        color = next(colors)

        if (str(aperture) != '100'):

            ax1.plot(values[1], values[2]/reference,
                     #color='red',
                     color=color,
                     #   marker='^',
                     #  linestyle='None',
                     label='ratio '+aperture+'/100')
            ax1.legend(loc='upper right')

        if (str(aperture) == '100'):
            pl.figure()
            title = str(img)+' ' + str(filters)
            print('title :', title, img)
            pl.title(title)
            pl.xlim([700, 1000])
            pl.ylabel('aper_flux ')
            pl.xlabel('wght (nm)')
            pl.plot(values[1], values[2],
                    #color='blue',
                    color=color,
                    label='flux ' + aperture)

            pl.legend(loc='upper right')
    pl.show()
