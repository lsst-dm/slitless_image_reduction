#!/usr/bin/env python 


import os, sys
import optparse
import astropy.io.fits as pyfits
from astropy import wcs
#import PythonPhot as pp
import numpy as np

import panstarrs as pan

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (13, 8) if False else (10, 6)

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()
 
    parser.add_option("-f","--file",default="../../data/c7675g0008o44.fits")
    parser.add_option("-x","--x0", type=float,default=1000.0)
    parser.add_option("-y","--y0", type=float,default=2000.0)
    parser.add_option("--xerr", type=float,default=10.0)
    parser.add_option("--yerr", type=float,default=10.0)
    parser.add_option("-r","--fwhm", type=float,default=100.0)

    parser.add_option("--doFake",  action="store_true", default=False)
    parser.add_option("-a","--famp", type=float,default=100.0) 

    opts, args = parser.parse_args()
  
    return opts

# Parse command line
opts = parse_commandline()

fitsfile = opts.file
name = fitsfile.split("/")[-1].replace(".fits","")

plotDir = './plots'
if not os.path.isdir(plotDir): os.mkdir(plotDir)
plotDir = './plots/%s'%name
if not os.path.isdir(plotDir):
    os.mkdir(plotDir)

ps1       = pan.panstarrs(fitsfile)
image    = ps1.OverscanSubtract_andTrim(gain=True)
hdr       = (ps1.fits[0].header).copy()
maskim = ps1.Mask() 

maskim = np.round(maskim)
image[maskim == 0] = 0.0
#image[maskim == 0] = np.nan
print np.nanmin(image), np.nanmedian(image), np.nanmax(image)

x0 = opts.x0
y0 = opts.y0
fwhm = opts.fwhm
xerr = opts.xerr
yerr = opts.yerr
famp = opts.famp

print 'fwhm', fwhm

x0index = int(np.round(x0))
y0index = int(np.round(y0))
fwhmindex = int(fwhm)
print 'fwhm', fwhmindex
xerrindex = int(np.round(xerr))
yerrindex = int(np.round(yerr))

if opts.doFake:
    pixels = np.arange(-fwhmindex*2,fwhmindex*2)
    [X,Y] = np.meshgrid(pixels,pixels)
    Z = famp*np.exp(-4*np.log(2) * (X**2 + Y**2) / fwhm**2)
    image[Y+y0index,X+x0index] = image[Y+y0index,X+x0index] + Z

fitsfile = os.path.join(plotDir,'image.fits')
hdu = pyfits.PrimaryHDU(image)
hdu.writeto(fitsfile,clobber=True)

vmin = np.nanmin(image)
vmax = np.nanmax(image)

ymax,xmax = image.shape

# Now let's see how we did
plt.figure()
plt.imshow(image,vmin=vmin,vmax=vmax,cmap='Greys_r', origin='lower')
plt.xlim([0,xmax])
plt.ylim([0,ymax])
# not bad!
# Save figure
plotName = os.path.join(plotDir,'image.png')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'image.eps')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'image.pdf')
plt.savefig(plotName)
plt.close()


sys.exit()


# first the sky background - let's just call the readnoise 0
skymod, skysig, skyskw = pp.mmm.mmm(image)

# This algorithm is not the best, helps to tune hmin a bit
hmin = 0.0
roundlim=[-10.0,10.0]
sharplim=[0.0,100.0]
xstar,ystar,flux,sharp,round =     pp.find.find(image,hmin,fwhm,roundlim=roundlim,sharplim=sharplim)

xstar,ystar = pp.cntrd.cntrd(image,xstar,ystar,fwhm,verbose=False)
index1 = np.where(np.abs(xstar-x0) < xerr)[0]
index2 = np.where(np.abs(ystar-y0) < yerr)[0]
index = np.intersect1d(index1,index2)
  
if len(index) == 0:
    print len(xstar), len(index)
    print "No stars found..."
    exit(0)

xstar = xstar[index]
ystar = ystar[index]
flux = flux[index]
sharp = sharp[index]
round = round[index]

a = np.vstack((np.around(xstar,1),np.around(ystar,1))).T
#a = np.vstack((xstar,ystar)).T
locs = np.vstack({tuple(row) for row in a})
xstar = locs[:,0]
ystar = locs[:,1]

gain = 1.0
mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = pp.aper.aper(image,xstar,ystar,phpadu=gain,apr=2*fwhm,zeropoint=0,skyrad=[3*fwhm,5*fwhm],exact=True)

vmin = np.nanmin(image)
vmax = np.nanmax(image)

# Now let's see how we did
plt.figure()
plt.imshow(image,vmin=vmin,vmax=vmax,cmap='Greys_r', origin='lower')
plt.plot(xstar,ystar,'o',ms=5,mfc='none',lw=2,mec='r')
plt.xlim([0,xmax])
plt.ylim([0,ymax])
# not bad!
# Save figure
plotName = os.path.join(plotDir,'image_sources.png')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'image_sources.eps')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'image_sources.pdf')
plt.savefig(plotName)
plt.close()

xmin = x0index - fwhmindex
xmax = x0index + fwhmindex
ymin = y0index - fwhmindex
ymax = y0index + fwhmindex
imagecut = image[ymin:ymax,xmin:xmax]

fitsfile = os.path.join(plotDir,'imagecut.fits')
hdu = pyfits.PrimaryHDU(imagecut)
hdu.writeto(fitsfile,clobber=True)

vmin = np.nanmin(imagecut)
vmax = np.nanmax(imagecut)
ymax,xmax = imagecut.shape

# Now let's see how we did
plt.figure()
plt.imshow(imagecut,vmin=vmin,vmax=vmax,cmap='Greys_r', origin='lower')
#plt.plot(xstar,ystar,'o',ms=5,mfc='none',lw=2,mec='r')
plt.xlim([0,xmax])
plt.ylim([0,ymax])
# not bad!
# Save figure
plotName = os.path.join(plotDir,'image_cut.png')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'image_cut.eps')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'image_cut.pdf')
plt.savefig(plotName)
plt.close()

xstar = xstar.flatten()
ystar = ystar.flatten()
mag = mag.flatten()
magerr = magerr.flatten()
flux = flux.flatten()
fluxerr = fluxerr.flatten()
sky = sky.flatten()
skyerr = skyerr.flatten()

magfile = os.path.join(plotDir,'magfile.dat')
fid = open(magfile,'w')
if isinstance(xstar, (list, tuple, np.ndarray)):
    for x,y,m,merr,f,ferr,s,serr in zip(xstar,ystar,mag,magerr,flux,fluxerr,sky,skyerr):
        fid.write('%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n'%(x,y,m,merr,f,ferr,s,serr))
else:
    fid.write('%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n'%(xstar,ystar,mag,magerr,flux,fluxerr,sky,skyerr))
fid.close()

