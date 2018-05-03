#!/usr/bin/env python 
'''
return plots of equivalent width

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys, re
import numpy as np
import pylab as pl
import toolbox as tb
import croaks
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import pandas as pd

def func(x, a, alpha):
    return a * x**alpha


objects =['HEN2-5',
          'HD185975',
                    'HD60753',
                    'HD14943',
                    'MUCOL',
                    'HD205905',
                    'HR7950',
                    'HD108344',
                    'HR9087']

filters=['NONE+RONCHI200',
         'RONCHI200+SEMROCK',
         'RONCHI200+z',
         'NONE+RONCHI400']

#objects=['HD14943']
#filters=['NONE+RONCHI400']
    
if __name__ == "__main__":

    airmass = []
    jd      = []
    wgth    = []
    ew      = []
    sew     = []
    img     = []
    seeing  = []
    zenith  = []
    ew_file = sys.argv[1]
    files   = sys.argv[2:]

   

    for obj in objects:
        for filt in filters:
            colors  = iter(cm.rainbow(np.linspace(0, 2, 200)))
            print obj, filt
            fig, axarr = pl.subplots(2, sharex=True)
            nb = 0
            for file in files:
                split = re.split('[_.]', file)
                Object = split[1]
                Filter = split[2]
                img    = re.split('/', split[0])[-1]
                if ((obj==Object) and (filt==Filter)):
                    nb+=1
                    #value, keys = tb.readlist(file)
                    values  = tb.readtxt(file, [
                        'pixel','aper_flux','center', 'psf_flux','sigma'])
                    pl.ylabel('aper_flux / psf_flux')
                    pl.xlabel('pixel')
                    title = str(obj)+' '+ str(filt)
                    print 'title :', title, img
                    pl.title(title)
                    color=next(colors)
                 
                    axarr[0].plot(values[0], values[1],
                            color=color,
                            linewidth=1,
                            #linestyle='None',
                            drawstyle='steps',
                            #capsize = 0,
                            label= str(img))
                    axarr[1].plot(values[0], values[1]/values[3],
                            color=color,
                            marker='o',
                            linestyle='None',
                            label= str(img))
                    pl.ylim(0.5,1.5)
             

            if(nb != 0):
                axarr[0].legend(loc='upper right')
                pl.show()
            pl.close('all')    
            pl.clf()
    sys.exit()
                
