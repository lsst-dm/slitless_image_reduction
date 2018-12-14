#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import re
import numpy as np
import pylab as pl
import exceptions
import __builtin__
import optparse
import scipy.optimize as optimization
import croaks
# ===========================================================


class fit_gains(Exception):
    pass
# ===========================================================


def sqr(a):
    return a*a


class fit_gains(object):
    def __init__(self, xmin, xmax):
        self.xmin = float(xmin)*2.0
        self.xmax = float(xmax)*2.0

    def func(self, x, a, b, c):
        return a*x*x + b*x + c

    ### Get variables from their names, no RO noise subtracted to the variance, the renorm is taken into account
    def get_mean_var(self, inputfile, ampli):
        data = croaks.NTuple.fromfile(inputfile)
        meanflatC = []
        varianceC = []
        reject_pixC = []
        residuC = []
        residu = (data[:]['residu'])
        reject_pix = (data[:]['reject_pix'])
        ampl = (data[:]['ampl'])
        meanflat = ((data[:]['mean1'])-(data[:]['ped1'])) + ((data[:]['mean2'])-(data[:]['ped2']))
        variance = (data[:]['var']) / ((data[:]['alpha']))
        for a, b, c, d, e in zip(meanflat, variance, ampl, reject_pix, residu):
            if (int(ampli) == int(c)):
                if (a > (self.xmin)):
                    meanflatC.append(float(a))
                    varianceC.append(float(b))
                    reject_pixC.append(float(d))
                    residuC.append(float(e))
        return meanflatC, varianceC, reject_pixC, residuC

    def fitting_interval(self, meanflat, variance, reject_pix, residu, ampli):
        c = zip(meanflat, variance, reject_pix, residu)
        dtype = [('mean', float), ('var', float), ('rej', float), ('res', float)]
        a = np.array(c, dtype=dtype)
        aa = np.sort(a, order='mean')
        frac_of_ptc = 0.0

        if (self.xmax == 0.0):
            tronc = 3 # maximum - $tronc points (usual safe value = 3)
            for item, i in enumerate(aa[:len(aa)-1]):
                if(i[1] > aa[item+1][1]):
                    self.xmax = aa[item-tronc][0]
                    frac_of_ptc = (aa[item-tronc][0]) / aa[item][0]
                    break
        del meanflat[:]
        del variance[:]
        del reject_pix[:]
        del residu[:]

        print('self.xmax', self.xmax)
        for i, j, k, l in aa:
            if (i <= self.xmax):
                meanflat.append(i)
                variance.append(j)
                reject_pix.append(k)
                residu.append(l)
        return meanflat, variance, reject_pix, residu

    def determine_gain(self, meanflat, variance, reject_pix, residu, ampli):
        GAIN_HEADER = 5.0
        #nb_pixel    = 1E6 # approx a 5% je pense
        nb_pixel = np.ones((len(variance)))*920400.0
        nb_pixel = nb_pixel-reject_pix
        x0 = np.array([-0.00001, 0.0001, 0.0001])
        zeros = np.zeros((len(variance)))
        #incertitude de la variance :
        meanflat = np.array(meanflat)
        variance = np.array(variance)
        sigma = 2 * np.array(variance) / np.sqrt(nb_pixel)
        param = 3
        fitR = optimization.curve_fit(self.func, meanflat, variance, x0, sigma)
        fitGain = fitR[0]
        covar = fitR[1]
        a = fitGain[0]
        b = fitGain[1]
        c = fitGain[2]
        c_error = np.sqrt(covar[2][2])
        b_error = np.sqrt(covar[1][1])
        a_error = np.sqrt(covar[0][0])
        khi2 = sqr((variance - (fitGain[0]*meanflat*meanflat + fitGain[1]*meanflat + fitGain[2]))/sigma)
        khi2 = khi2.sum()
        KN = khi2 / ((len(sigma)) - param)
        g = 1/fitGain[1]
        sg = b_error/(fitGain[1]*fitGain[1])
        #print "Gain = ", g , " +/- ", sg
        #print "sigma = ", sigma.mean(), " variance = ", variance.mean()
        fig = pl.figure(1)
        p1 = pl.errorbar(meanflat, variance, yerr=sigma, xerr=zeros, fmt='b^', markersize=6, label="Mes. sx")
        #fitquad = pl.polyval(a*meanflat*meanflat + b*meanflat + c, meanflat)
        #fitquad = pl.polyval(2*meanflat, meanflat)
        meanflat = np.sort(meanflat)
        ratio = (((g*2E4) - (a*(2E4*g)*(2E4*g) + b*(2E4*g) + c)) / (a*(2E4*g)*(2E4*g) + b*(2E4*g) + c))
        print("gain = ", g, " khi2/NdF = ", KN, " ratio = ", ratio)
        if (ampli <= 7): ## various checks
            pl.plot(meanflat, a*meanflat*meanflat + b*meanflat + c, 'r')
            #pl.plot(meanflat,b*meanflat + c,'b')
        #pl.plot(meanflat,( meanflat/b + (1/(1+meanflat*c))), 'r')
            #pl.plot(meanflat/2.0,(variance-(a*meanflat*meanflat + b*meanflat + c))/(a*meanflat*meanflat + b*meanflat + c),'r')
        pl.show()
        out = [float(ampli), GAIN_HEADER, g, sg, c, c_error, a, a_error, KN, self.xmax, ratio]
        return out

    def dump_gain_and_parameters(self, inputfile, ampli):
        meanflat, variance, reject_pix, residu = self.get_mean_var(inputfile, ampli) # Get the data
        meanflat, variance, reject_pix, residu = self.fitting_interval(
            meanflat, variance, reject_pix, residu, ampli) # select the span of the fit
        return self.determine_gain(meanflat, variance, reject_pix, residu, ampli)
