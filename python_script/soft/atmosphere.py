#!/usr/bin/env python 
'''
return wght calibrated spectrum + forward model + EW measurements

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys, re
import numpy as np
import pylab as pl
import toolbox as tb
import matplotlib.cm as cm
from scipy.optimize import curve_fit
import scipy.interpolate as interp
import sed as Sed
import lines as Lines
import argparse
import logging


'''-------------------------------------'''

class Modeling(object):
    def __init__(self, target, ronchi,
                 convolve = 2., 
                 directory = './',
                 plot      = False,
                 pix2wgth = 2.1,
                 keys = '',
                 telescope = 'open',
                 seeing_at_wght = [2.3e-06,  -3.9e-03, 3.6],
                 **kwargs):
        self.target               = target
        self.keys                 = keys 
        self.ronchi               = ronchi
        self.directory            = directory
        self.plot                 = plot
        self.pix2wgth_coef        = pix2wgth 
        self.convolve             = convolve
        self.telescope            = telescope
        self.Init(seeing_at_wght)


    def Init(self, seeing_at_wght):
        self.sed = Sed.SED(self.target, self.convolve,
                           pix2wgth       = self.pix2wgth_coef,
                           ronchi         = self.ronchi,
                           plot           = self.plot,
                           seeing_at_wght = seeing_at_wght)
   
        self.sed_wavelength, self.sed_flux = self.sed.smoothedSed()
      

        
    '''Tatmo = Sobs /(SED*Ttel)'''
    def Tatmo(self, wavelength, signal, mkcalib=True, name = '', flux='', **kwargs):
        w, Ttel = tb.telescope_T(self.telescope)
        Ttel    = interp.griddata(w, Ttel, self.sed_wavelength)
        SedTel  = self.sed_flux * Ttel
        norm    = np.nanmax(SedTel)
        print norm
        SedTel  = np.array(SedTel)/norm
        signal  = np.array(signal)           
        signal  = interp.griddata(wavelength, signal, self.sed_wavelength)
        norm    = np.nanmax(signal)
        signal  = signal/norm
        tatmo   = signal / SedTel
        w, tatmo= tb.clean(self.sed_wavelength, tatmo)
        #for i, j in zip(w, tatmo):
        #   print i, j

        tb.DumpFile(self.keys,
                    ['w','tatmo'],
                    [w, tatmo],
                    os.path.join(self.directory,str('Tatmo'+name+'.list')))
        print 'writing : ', os.path.join(self.directory,str('Tatmo'+name+'.list')) 
        if (self.plot==True):
            fig2 = pl.figure()
            pl.plot(self.sed_wavelength, Ttel, label='Telescope Throughput')
            pl.xlabel('wavelength (nm)')
            pl.ylabel('Telescope Throughput')
            pl.legend()
            fig3 = pl.figure()
            pl.plot(self.sed_wavelength, SedTel, label='SED * Ttel')
            pl.plot(self.sed_wavelength, signal, label='obs')
            pl.xlabel('wavelength (nm)')
            pl.ylabel('SED * Ttel')
            pl.legend()
            fig4 = pl.figure()
            pl.plot(w, tatmo, label='Tatmo from ' +str(flux))
            pl.xlabel('wavelength (nm)')
            pl.ylabel('Sobs /(SED*Ttel)')
            pl.legend()
            pl.show()
      

'''-------------------------------------'''



'''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
'''to dispersion direction (return poly2 fit parameters)                         '''
def SeeingAtWght(w, psf_sigma):
    clean          = np.array([w, psf_sigma]).transpose()
    clean          = clean[(clean[:,0]>400) & (clean[:,0]<900)]
    seeing_at_wght = np.polyfit(clean[:,0], clean[:,1], deg=2)
    return seeing_at_wght 



def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "process slitless images"
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="analyze spectrum")
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
    parser.add_argument('-e',"--flux", 
		        help = "flux estimator, default is aper", 
		        default='aper_flux')
    parser.add_argument('-v',"--verbose", 
		        help = "verbose", 
		        action='store_true')
    parser.add_argument('-f',"--file", type=str, 
	                help = "name of spectrum file", 
	                default='rawspectrum.list')
    parser.add_argument('-n',"--name", type=str, 
	                help = "output files identifiers", 
	                default='+1')
    parser.add_argument('-i',"--imgs", nargs='+', type=str, 
	                help = "list of spectrum to be analyzed", 
	                default=None)
    args = parser.parse_args()
    return args




if __name__ == "__main__":
   
    # Run the atmo after the recalibration of the disperison relation.

    # Dump a spectrum.list with the old and new dispersion, and a model column.
    # Also dump a pdf comparing the template with the recalibrated observation.
    # From the recalib. dump new EW estimate from fix edges positions.

    args       = grabargs()
    plot       = args.plot
    verbose    = args.verbose
    input_reps = args.imgs
    file       = args.file
    name       = args.name
    flux       = args.flux
    
    if(verbose):
        Level=logging.getLogger().setLevel(logging.DEBUG)
        logging.debug('DEBUG mode')
    else:
        Level=logging.getLogger().setLevel(logging.INFO)
        logging.info('INFO mode')
     
  
    tb.ProdInstall()

    template_path = os.environ['TEMPLATE_PREFIX']
    if not os.path.exists(template_path):
        os.makedirs(template_path)
  
    for input_rep in input_reps :
        logging.info('\n' + 'Analyzing : '+ input_rep)

        
        '''observed spectrum'''
        input_file  =  os.path.join(input_rep, file) 
        if not os.path.isfile(input_file):
            continue

        print 'reading ', input_file
        dict, array, names = tb.readcat(input_file)
        target  = ' '.join(dict.get('OBJECT'))
        print target
        if not target:
            continue

        filters    = ' '.join(dict.get('FILTERS'))
        if (str(filters).find('RG715') >= 0):
            print 'RG715 filter in place'
            tel_t = 'rg715'
        else:
            print 'Telescope set up open'
            tel_t = 'open'
            
        seeing      = float(dict.get('SEEING')[0]) # This is currently not a good measurement
        pix2wgth    = float(dict.get('PIX2WGTH')[0])
        ronchi      = dict.get('RONCHI')  [0]
        array           = np.array(array.tolist())#np.array(array).transpose()

        array           = array[(array[:,1]<=1100)] # Need a better model for the fitting

        pixel           = array[:,0]
        w               = array[:,1]
        aper_flux       = array[:,2]
        psf_gauss_flux  = array[:,3]
        psf_gauss_sigma = array[:,4]
        psf_gauss_mu    = array[:,5]
        psf_moffat_flux = array[:,6]
        psf_moffat_x0   = array[:,7]
        psf_moffat_gamma = array[:,8]
        integralGM      = array[:,9]
        integralG       = array[:,10]
        amplitude_0     = array[:,11]
        x_0_0           = array[:,12]
        gamma_0         = array[:,13]
        alpha_0         = array[:,14]
        amplitude_1     = array[:,15]
        mean_1          = array[:,16]
        stddev_1        = array[:,17]

        
      
        

        '''Select the best flux estimator'''
        #best_flux = psf_gauss_flux #psf_voigt_flux
        #best_flux = flux #aper_flux #psf_voigt_flux
        if flux == 'aper_flux':
            best_flux = aper_flux
         
        if flux == 'psf_gauss_flux':
            best_flux = psf_gauss_flux
            
        if flux == 'psf_moffat_flux':
            best_flux = psf_moffat_flux


        if flux == 'psf_GM_flux':
            best_flux = integral_GM

        if flux == 'flux_deco':
            best_flux = array[:,15]

        print 'flux :', flux
        if plot:
            pl.plot(w, best_flux)
            pl.show()
            
        '''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
        '''to dispersion direction (return poly2 fit parameters)                         '''
        seeing_at_wght = SeeingAtWght(w, psf_gauss_sigma)

        
        
        ''' Load CalSpec reference spectrum based on target name   '''
        ''' Accounting for atmospheric transmission and resolution '''
        ''' Atmospheric transmission from forward modeling         '''
     
        convolve = pix2wgth * seeing # in nm # constant approximation, deprecated
        
        mod = Modeling(target, ronchi,
                       convolve  = convolve,
                       directory = input_rep,
                       plot      = plot,
                       keys      = dict,
                       pix2wgth  = pix2wgth,
                       telescope = tel_t,
                       seeing_at_wght = seeing_at_wght)

       
        #name  = re.sub('/','_', input_rep[:-2])
        name  = '_'.join((input_rep.split('/'))[-3:-1])
        mod.Tatmo(w, best_flux, mkcalib=True, name = name, flux = flux)
