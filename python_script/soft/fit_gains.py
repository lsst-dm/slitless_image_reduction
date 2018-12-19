#!/usr/bin/env python
from __future__ import print_function
from builtins import str
from builtins import next
import os
import sys
import re
import numpy as np
import pylab as pl
import scipy.optimize as optimization
import croaks
import toolbox as tb
import argparse
import matplotlib.cm as cm


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "quadratic of gain on ptc"
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract spectrum from raw image")

    parser.add_argument('-i', "--img", type=str,
                        help="object name if not in fits header",
                        default=None)

    args = parser.parse_args()
    return args


def func(x, a, b, c):
    return a*x*x + b*x + c


def determine_gain(meanflat, variance, ampli, color):
    GAIN_HEADER = 5.0
    nb_pixel = 1E6 # approx a 5% je pense
    x0 = np.array([-0.00001, 0.0001, 0.0001])
    zeros = np.zeros((len(variance)))
    #incertitude de la variance :
    meanflat = np.array(meanflat)
    variance = np.array(variance)
    sigma = 2 * np.array(variance) / np.sqrt(nb_pixel)
    param = 3
    fitR = optimization.curve_fit(func, meanflat, variance, x0, sigma)
    fitGain = fitR[0]
    covar = fitR[1]
    a = fitGain[0]
    b = fitGain[1]
    c = fitGain[2]
    c_error = np.sqrt(covar[2][2])
    b_error = np.sqrt(covar[1][1])
    a_error = np.sqrt(covar[0][0])
    khi2 = ((variance - (fitGain[0]*meanflat*meanflat + fitGain[1]*meanflat + fitGain[2]))/sigma)**2
    khi2 = khi2.sum()
    KN = khi2 / ((len(sigma)) - param)
    g = 1/fitGain[1]
    sg = b_error/(fitGain[1]*fitGain[1])
    print('Ampli : ', ampli, " gain = ", g, " +/- ", sg, ' e/ADU')
    #print "sigma = ", sigma.mean(), " variance = ", variance.mean()
    fig = pl.figure(1)
    p1 = pl.errorbar(meanflat, variance, yerr=sigma, xerr=zeros, fmt='^', markersize=6, color=next(colors))
    #fitquad = pl.polyval(a*meanflat*meanflat + b*meanflat + c, meanflat)
    #fitquad = pl.polyval(2*meanflat, meanflat)
    meanflat = np.sort(meanflat)
    #ratio = (((g*2E4) - (a*(2E4*g)*(2E4*g) + b*(2E4*g) + c)) / (a*(2E4*g)*(2E4*g) + b*(2E4*g) + c))
    print(" khi2/NdF = ", KN)
    pl.plot(meanflat, a*meanflat*meanflat + b*meanflat + c, color='k',
            label='amp ' + str(ampli)+' g = '+str(round(g, 2))+' (e/ADU)')
    #pl.plot(meanflat,b*meanflat + c,'b') # show linear part
    #pl.plot(meanflat,( meanflat/b + (1/(1+meanflat*c))), 'r')
    pl.xlabel('flux (ADU)')
    pl.ylabel('variance')
    pl.legend(loc='upper left')
    pl.title('Quadratc fit of the PTC')
    out = [float(ampli), GAIN_HEADER, g, sg, c, c_error, a, a_error, KN]
    return out


if __name__ == '__main__':
    args = grabargs()
    file = args.img
    flux_min = 0.0
    flux_max = 58000.0 * 2. # it's a pair!
    values = tb.readtxt(file, ['ampl', 'mean_sum', 'rms_diff'])
    values = np.array(values).transpose()
    amps = ('11', '12', '21', '22')
    colors = iter(cm.rainbow(np.linspace(0, 1, 4)))
    for amp in amps:
        ptc = values[(values[:, 0] == int(amp)) & (values[:, 1] >= flux_min) & (values[:, 1] <= flux_max)]
        meanflat = ptc[:, 1]
        variance = ptc[:, 2]**2
        determine_gain(meanflat, variance, amp, colors)

    pl.show()
