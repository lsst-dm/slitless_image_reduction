#!/usr/bin/env python 
'''
splines to extract continuum
Author : Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys, re
import numpy as np
import pylab as pl
import toolbox as tb
import scipy.interpolate as interp
import itertools
import smooth 
import argparse



'''fit using features position found in both the obs. and the template'''


'''Return the coefficients to adjust the wght solution on the template'''
'''Aslo return the list of features detected on both the template and the obs'''
def Fit(Lm, Lo, plot=False, **kwargs):
    degree = 2
    if (len(Lm) <= degree):
        degree = int(len(Lm)-1)
        print "Changed degree to : ", degree
        
    fit_coef    = np.polyfit(Lo, Lm, deg=degree)
    points      = np.linspace(100, 1500 ,100)
    dispersion  = np.polyval(fit_coef, Lo)
    dispersion2 = np.polyval(fit_coef, points)
    if plot is True :

        fig, ax1 = pl.subplots()
        pl.xlabel('obs (nm)')
        pl.ylabel('template (nm)')  
        ax1.plot(Lo, Lm, 'r^')
        ax1.plot(points, dispersion2, color = 'k')
        pl.show()
      
        fig, ax2 = pl.subplots()
        ax2.plot(Lo, Lm-dispersion, 'r^')
        pl.axhline(y=0., linewidth=1, color = 'k')
        pl.xlabel('wght (nm)')
        pl.ylabel('residuals')
        pl.show()
    print 'rms residuals = ' , np.std(Lm-dispersion)
    return fit_coef, Lm
    



'''Extraction of EW, fitting 2-d polynom on edges'''
def equWidth2(wavelength, flux, list, slide = 0., plot=False, **kwargs):
    e1a = list[1][0] + slide
    e1b = list[1][1] + slide
    e2a = list[2][0] + slide
    e2b = list[2][1] + slide
    data    = np.array([wavelength, flux]).transpose() 
    edge1   = data[(data[:,0]>= e1a) &  (data[:,0]<= e1b)] 
    edge2   = data[(data[:,0]>= e2a) &  (data[:,0]<= e2b)] 

    fit1 = np.polyfit(edge1[:,0], edge1[:,1], deg=2)
    fit2 = np.polyfit(edge2[:,0], edge2[:,1], deg=2)


    e1 = np.polyval(fit1, e1b) 
    e2 = np.polyval(fit2, e2a)
    
    fit = np.polyfit([e1b, e2a], [e1, e2], deg=1)
    wght      = np.linspace(e1b, e2a, int(e2a-e1b))
    continuum = np.polyval(fit, wght)

    signal    = interp.griddata(wavelength, flux, wght)

    spacing = wght[1]-wght[0]
    equ_w   = np.trapz(1- signal / continuum, dx=spacing)
                           
    if plot:
        pl.figure()
        pl.plot(wavelength, flux, label='signal')
        pl.plot(edge1[:,0], edge1[:,1], 'r^', label='edge1')
        pl.plot(edge2[:,0], edge2[:,1], 'r^', label='edge2')

        points1  = np.linspace(e1a , e1b, 20)
        out1     = np.polyval(fit1, points1)
        pl.plot(points1, out1, linewidth=1, color = 'k')
        points2  = np.linspace(e2a, e2b, 20)
        out2     = np.polyval(fit2, points2)
        pl.plot(points2, out2, linewidth=1, color = 'k')
        pl.plot(wght, continuum, linewidth=1, color = 'g')
        pl.legend()
        pl.show()
          
    segment = (e2a + e1b)/2.
    return segment, equ_w





def fromSimu(file):
    values = tb.readtuple(file)
    data = np.array(values).transpose()
    wght = [float(i) for i in data[0]]
    flux = [float(i) for i in data[1]]
    data = np.array([wght, flux]).transpose()
    try:
        airmass = float((((os.path.basename(file)).split('_'))[6]).strip('z'))/1000.
    except:
        airmass = 0

    try:
        pwv = float((((os.path.basename(file)).split('_'))[7]).strip('wv'))/100.
    except:
        pwv = 0
    Object  = 'ref'
    jd      = 0.
    parallactic_angle = 0.
    return data, airmass, Object, jd, parallactic_angle, pwv

def fromObs(file):
    dict, values, names = tb.readcat(file)
    meanseeing = meanSeeing(values.field('w'), values.field('psf_gauss_sigma'))
    #print 'Dict : ' ,dict.items()
    latitude = ' '.join(dict.get('LATITUDE'))
    ha       = ' '.join(dict.get('HA'))
    dec      = dict.get('DEC')
    parallactic_angle = tb.ParallacticAngle(latitude, ha, dec)
    parallactic_angle = parallactic_angle[0]
    print 'parallactic_angle ', parallactic_angle
    Object   = ' '.join(dict.get('OBJECT'))
    jd       = dict.get('JD')[0]
    airmass  = dict.get('AIRMASS')[0]
    #data = np.array([values.field('w'), values.field('aper_flux')]).transpose()
    data = np.array([values.field('w'), values.field('psf_gauss_flux')]).transpose()
    print 'flux estimator is Gauss PSF'
    return data, airmass, Object, jd, parallactic_angle, meanseeing


'''
HD14943 absoption lines (nm):
1 410.3
2 434.2
3 486.2
4 656.5
5 762.0
'''

EW =[['1' , [401,406], [414,419], 410.3 ] ,
     ['2' , [425,430], [437,442], 434.2 ] ,
     ['3' , [470,480], [492,502], 486.2 ] ,
     ['4' , [640,650], [665,675], 656.5 ] ,
     ['5' , [735,745], [780,785], 762.0 ] ]


def meanSeeing(wght, sigma):
    data = np.array([wght, sigma]).transpose()
    seeing = data[:,1][(data[:,0]>= EW[1][1][1]) & (data[:,0]<= EW[1][2][1])]
    return np.mean(seeing)
    


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum and extract EW"
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="analyze spectrum")
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
    parser.add_argument('-f',"--outfile", type=str, 
	                help = "name of output file", 
	                default='regression.list')
    parser.add_argument('-s',"--seeing", type=str, 
	                help = "seeing from file", 
	                default=None)
    parser.add_argument('-t',"--Type", type=str, 
	                help = "either obs or simu", 
	                default=None)
    parser.add_argument('-i',"--files", nargs='+', type=str, 
	                help = "input files", 
	                default=None)
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args     = grabargs()
    plot     = args.plot
    outfile  = args.outfile
    inputype = args.Type #''' Either obs or simu '''
    files    = args.files 
    seeing   = args.seeing


    if ((inputype=='simu') and (seeing is not None)):
        dict, values, names = tb.readcat(seeing)
        seeing_x = values.field('w')
        seeing_y = values.field('psf_gauss_sigma')
        pix2wght = float(dict.get('PIX2WGTH')[0])
        mean_seeing = meanSeeing(seeing_x, seeing_y)
        print 'MEAN SEEING ' , mean_seeing
    

    outlist = []
    for file in files:
        print 'opening ', file
        if (inputype=='obs'):
            data, airmass, Object, jd, parallactic_angle, mean_seeing = fromObs(file)       
            pwv = 0.
        if (inputype=='simu'):
            data, airmass, Object, jd, parallactic_angle, pwv = fromSimu(file)


        wmin = 360.
        wmax = 1000.
        data = data[(data[:,0]>= wmin) & (data[:,0]<= wmax)]
   
        ''' 
        Smoothing the atmospheric curve to match the resolution of the observation
        '''
        if ((inputype=='simu') and (seeing is not None)):
             sm   = smooth.smooth(data[:,0], data[:,1], seeing_x, seeing_y,
                                      pix2wgth = pix2wght,
                                      plot = plot)
             data[:,0], data[:,1] = sm.smoothCurve()
        if seeing is None and inputype is 'simu':
            mean_seeing = 0.
   

    
        '''
        --> Fitting a polynom on edges 
        --> measurement of the EW
        --> +/- 10 nm sliding
        --> picking the max
        --> fitting a parabola on data-continuum and retrieving the minimum
        '''
        for i in EW:
            minimize = []
            for slide in range(-10,11,1):
            
                print 'EW : ', i, ' slide = ', slide
                ew_w, ew_f = equWidth2(data[:,0],data[:,1], i,
                                       slide=slide,
                                       plot=False)         
                print  ew_w, ew_f
                minimize.append([slide, ew_w, ew_f])


            search = list(map(list, zip(*minimize)))
            MAX = max(search[2])
            for s in minimize:
                print s
                if s[2]==MAX:
                    f_slide = s[0]
                    center  = s[1]

            print f_slide, center

            ew_w, ew_f = equWidth2(data[:,0],data[:,1], i,
                                   slide = f_slide,
                                   plot = plot)         
            
            print ew_w, ew_f 
            outlist.append([Object, i[3], jd, airmass,parallactic_angle, ew_w, ew_f, mean_seeing, pwv])


        names =['object', 'expected', 'jd', 'airmass', 'parallactic','ew_w','ew_f', 'mean_seeing', 'pwv']
        tb.DumpTuple(names,
                     zip(*outlist),
                     outfile)
        print 'writing : ', outfile
     


        ''' Refit the dispersion relation by matching absorption features '''
        ''' between the observation and a template                '''
        search = list(map(list, zip(*outlist)))
        Lm = search[1]
        Lo = search[5]
        fit_coef, ew_position = Fit(Lm, Lo, plot=plot)
      
        ''' Shifting the initial wght calibration ''' 
        wght_recal = np.polyval(fit_coef, data[:,0])

        if plot :
            pl.plot(data[:,0], data[:,1], label='initial')
            pl.plot(wght_recal, data[:,1], label='wght recalibrated')
            pl.xlabel('wght (nm)')
            pl.ylabel('best_flux (ADU)')
            pl.legend()
            pl.show()

        


        '''writing the wght_recal to file'''
        dict, values, names = tb.readcat(file)
        wght_recal = np.polyval(fit_coef, values.field('w'))
        val = np.array(values)
        from numpy.lib.recfunctions import append_fields
        z = append_fields(val, 'y', wght_recal)
        outcat = 'test.list'
        print 'writing ', outcat
        names = names + ['wght_recal']
        print 'output entries : ', names
        tb.DumpTuple(names, zip(*z), outcat)

    
