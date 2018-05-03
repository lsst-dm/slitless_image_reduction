#!/usr/bin/env python 
'''
fit a distribution, either Rayleigh or double gauss

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys, re
import numpy as np
import pylab as pl
import toolbox as tb

from scipy.optimize import curve_fit
from astropy.modeling import models, fitting
import argparse

# Define model function to be used to fit to the data above:
def rayleigh(x, *p):
    A, sigma = p
    return A*x*np.exp(-(x)**2/(2.*sigma**2)) / sigma**2



def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum and extract EW"
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="fit a distribution")
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
    parser.add_argument('-d',"--distrib", type=str, 
		        help = "choose distribution to fit", 
		        default='R')
    parser.add_argument('-m',"--min", type=float, 
	                help = "interval minimum", 
	                default=0.)
    parser.add_argument('-a',"--max", type=float, 
	                help = "interval maximum", 
	                default=1000.)
    parser.add_argument('-v',"--value", type=str, 
	                help = "name of observable", 
	                default='TQV')
    parser.add_argument('-s',"--site", type=str, 
	                help = "name of site", 
	                default='CTIO')
    parser.add_argument('-i',"--file", type=str, 
	                help = "input file", 
	                default=None)
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args    = grabargs()
    file1   = args.file
    plot    = args.plot
    name    = args.value
    site    = args.site
    method  = args.distrib
    MIN     = args.min
    MAX     = args.max
    print MIN, MAX
    [hour, jd, value] = tb.readtxt(file1, ['entry', 'jd', 'value'])
    base = 24. # if sampling is hourly
    
    if plot:
        fig = pl.figure(1)
        pl.plot(hour / base, value, label=site+' 2017')
        pl.xlabel('time (days)')
        pl.ylabel(name)
        pl.legend()
      

    
    value = value[(value<MAX)&(value>MIN)]

    fig = pl.figure(2)
    n, bins, patches = pl.hist(value, 50,
                               normed=True,
                               facecolor='g', alpha=0.75)

    bin_width = ( bins[1]- bins[0])
    center = bin_width/2.
    x = bins[0:-1]+center

    print 'histogram integral = ', sum(n)*bin_width
    
    if method=='R':
        amplitude  = 1.
        sigma = 2.
        p0 = [amplitude, sigma]   
        coeff, var_matrix = curve_fit(rayleigh, x, n,
                                      p0=p0)
        amplitude = coeff[0]
        sigma     = coeff[1]
        print "A*x*np.exp(-(x)**2/(2.*sigma**2)) / sigma**2"
        print "A, sigma = ", coeff
        fit = rayleigh(x, *coeff)
        distrib = "Rayleigh"
        pl.text(10, 0.1,
                'A, sigma =  %s'%(coeff),
                color='red', fontsize=8)
   
                
    if method=='G':
        gg_init = models.Gaussian1D(amplitude=1, mean=260, stddev=10.) + models.Gaussian1D(amplitude=1, mean=290, stddev=10.)
        fitter = fitting.SLSQPLSQFitter()
        gg_fit = fitter(gg_init, x, n)
        print gg_fit
        fit = gg_fit(x)
        distrib = "Two gaussian"
        pl.text(220, 0.02,
                '%s \n %s'%(gg_fit.param_names, gg_fit.parameters),
                color='red', fontsize=8)
   
   
    pl.plot(x, fit, 'b-', label=distrib+' distribution')
    pl.xlabel(name)
    pl.legend()
    pl.show()


    outname = name+'_distrib'
    names   = ['bin', 'value']
    tb.DumpTuple(names,
                [x,nTO3_distrib.list*bin_width],
                 outname+'.list')
    print 'writing : ', outname+'.list'
     
