#!/usr/bin/env python
'''
analyze ptc

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''


import os
import sys
import astropy.io.fits as pf
import glob
import telinst as instru
import pylab as pl
import pylab
import argparse
import croaks
import numpy as np
from astropy.stats import sigma_clip


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "process slitless images"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract spectrum from raw image")

    parser.add_argument('-i', "--img", type=str,
                        help="object name if not in fits header",
                        default=None)

    args = parser.parse_args()
    return args


def readtxt(inputfile, dic=[]):
    data = croaks.NTuple.fromtxt(inputfile)

    value = []
    for name in dic:
        value.append(np.array(data[:][name]))
    return value, data.keys


# correlations by  fft :

def compute_corr_fft(diff, w, maxrange, fft_shape):
    # check that the zero padding implied by "fft_shape"
    # is large enough for the required correlation range
    assert(fft_shape[0] > diff.shape[0]+maxrange+1)
    assert(fft_shape[1] > diff.shape[1]+maxrange+1)
    tim = np.fft.rfft2(diff*w, fft_shape)
    tmask = np.fft.rfft2(w, fft_shape)
    pcov = np.fft.irfft2(tim*tim.conjugate())
    pmean = np.fft.irfft2(tim*tmask.conjugate())
    pcount = np.fft.irfft2(tmask*tmask.conjugate())
    l = []
    # (dy,dx) = (0,0) has to be first
    for dy in range(maxrange+1)+range(-1, -maxrange-1, -1):
        for dx in range(0, maxrange+1):
            npix = pcount[dy, dx]
            cov = pcov[dy, dx]/npix-pmean[dy, dx]*pmean[-dy, -dx]/(npix*npix)
            if (dx == 0 and dy == 0):
                var = cov
            l.append("%d %d %f %f %d"%(dx, dy, var, cov, npix))
    return l


def Clip(img):
    mask = sigma_clip(img, sigma=4, iters=5, cenfunc=np.mean, copy=True)
    return mask


def Mask(inst, img, amp):
    img = inst.IlluRegion(img, amp)
    return Clip(img)


if __name__ == '__main__':
    args = grabargs()
    ptc = args.img
    [img1, img2], keys = readtxt(ptc, ['img1', 'img2'])

    outfile = 'ptc_var.list'
    out = open(outfile, 'w')
    out.write("""# ampl : 
# exptime :
# mean_sum : 
# rms_sum :
# mean_diff :
# rms_diff : 
# end 
""")

    number = 1
    for refname1, refname2 in zip(img1, img2):
        print refname1, refname2
        inst = instru.telinst(refname1)
        exptime = inst.header.get('EXPTIME')
        data1 = inst.Image(refname1)
        data2 = inst.Image(refname2)
        outimg1 = inst.OverscanSubtract_andTrim(data1)
        outimg2 = inst.OverscanSubtract_andTrim(data2)
        diff = outimg1 - outimg2
        hdr = (inst.header).copy()
        info = "Image is a stack of :" + refname1 + refname2
        hdr.add_comment(info)
        outim = 'img'+str(number)+'.fits'
        print 'outimg ', outim, ' exptime = ', exptime
        pf.writeto(outim, diff, hdr, overwrite=True)
        number += 1
        for amp in inst.amp:
            '''Construct a mask'''
            mask1 = Mask(inst, outimg1, amp)
            mask2 = Mask(inst, outimg2, amp)
            masked_diff = Clip(mask1 - mask2)
            masked_sum = Clip(mask1 + mask2)
            diff_mean = np.mean(masked_diff)
            diff_rms = np.std(masked_diff)
            print 'amp, diff_mean diff_std = ', amp, diff_mean, diff_rms
            sum_mean = np.mean(masked_sum)
            sum_rms = np.std(masked_sum)
            out.write('%i %i %f %f %f %f\n' %
                      (int(amp), int(exptime), sum_mean, sum_rms, diff_mean, diff_rms))

    out.close()

    #print mask1.fill_value
    #pl.figure()
    #pl.subplot(211)
    #pl.imshow(masked_diff, cmap='gray')
    #pl.colorbar()
    #pl.subplot(212)
    #pl.imshow(img, cmap='gray')
    #pl.colorbar()
    #pl.show()
