#!/usr/bin/env python
from __future__ import print_function
import os
import numpy as np
import poloka.core as pc
import poloka.bfstuff as bf
import pyfits
# ===========================================================


class pixel_utils(Exception):
    pass
# ===========================================================


def append_header_fromfile(header):
    """
    Adds header from file
    """
    defaultheader = "header.txt"
    headerhdu = pyfits.Header.fromtextfile(defaultheader)
    header.update(headerhdu)


def change_datasec(extheader, extnum, X0, Y0, NX, NY):
    """
    Corrects extension DETSEC, DATASEC
    """
    parstringlow = '1:2002'
    parstringhigh = '4004:2003'
    colwidth = 512
    extheader['DATASEC'] = '[%s:%s,%s:%s]'%(X0, NX, Y0, NY)
    if extnum < 9:
        pdet = parstringlow
        si = colwidth*(extnum)
        sf = colwidth*(extnum-1)+1
    else:
        pdet = parstringhigh
        si = colwidth*(16-extnum)+1
        sf = colwidth*(16-extnum+1)
    extheader['DETSEC'] = '[{}:{},{}]'.format(si, sf, pdet)


class Rxx_flat(object):
    def __init__(self, nsig, bordersize, minpixcount, biassec='[531:536,29:1800]', datasec='[37:510,19:1970]', device="LSST", exptime_name="0", **kwargs):
        # mask parameters
        #(1) flag above a nsig threshold
        #(2) size of patch around a flagged pixel
        #(3) nb of flag in the patch to mask the pixel
        self.device = device
        self.nsig = nsig
        self.bordersize = bordersize
        self. minpixcount = minpixcount
        self.bias = []
        self.overscan_mean_and_var = []
        self.BIASSEC = biassec
        self.DATASEC = datasec
        self.exptime_name = exptime_name
        self.reftime = self.exptime_name

    #(1) It changes keywords that are used by poloka.core routines
    def correct_dataAndBiasSec(self, extheader, extnum):
        """
        Corrects extension DETSEC, DATASEC, BIASSEC and DETSIZE
        """
        parstringlow = '1:2002'
        parstringhigh = '4004:2003'
        colwidth = 512
        extheader['DETSIZE'] = '[1:4096,1:4004]'
        extheader['DATASEC'] = self.DATASEC
        extheader['BIASSEC'] = self.BIASSEC
        if extnum < 9:
            pdet = parstringlow
            si = colwidth*(extnum)
            sf = colwidth*(extnum-1)+1
        else:
            pdet = parstringhigh
            si = colwidth*(16-extnum)+1
            sf = colwidth*(16-extnum+1)

        if(self.device == 'ITL'):
            extheader['DETSEC'] = '[4:512,1:2000]'
            extheader['DETSIZE'] = '[1:4336,1:4044]'
        else:
            extheader['DETSEC'] = '[{}:{},{}]'.format(si, sf, pdet)
        print('DETSEC ', extheader['DETSEC'])
        print('DATASEC ', extheader['DATASEC'])
        print('BIASSEC ', extheader['BIASSEC'])
        print('DETSIZE ', extheader['DETSIZE'])

    def correct_Sec_Subaru(self, extheader):
        extheader['DETSIZE'] = '[1:2048,1:4177]'
        extheader['DATASEC'] = self.DATASEC
        extheader['BIASSEC'] = self.BIASSEC
        illu = self.DATASEC
        b = "'["
        for char in b:
            illu = illu.replace(char, "")
        illu = illu.replace(':', ' ').replace(',', ' ').split()
        namp = 4
        dist = []
        for i in range(0, namp):
            dist.append(np.absolute(512*i-int(illu[0])))
        smallest = dist.index(min(dist))
        Iamp = smallest + 1
        incr = 512 * (Iamp-1)
        Xmin = 1 + incr
        Ymin = 1
        Xmax = 512 + incr
        Ymax = 4177
        extheader['DETSEC'] = '[%s:%s,%s:%s]'%(Xmin, Xmax, Ymin, Ymax)

    #(1) It strips area keywords
    def getsec(self, val):
        field = map(str, val.strip('[]').split(','))
        i_in = int((field[0].split(':'))[0])
        i_fin = int((field[0].split(':'))[1])
        j_in = int((field[1].split(':'))[0])
        j_fin = int((field[1].split(':'))[1])
        bs = [i_in, i_fin, j_in, j_fin]
        return bs

    #(1) it return time in sec since 1.1.1970 minus 1407000000 sec.

    def get_reftime(self, hdr):
        if (self.exptime_name == "0"):
            if "CTIME" in hdr:
                ref = int(hdr['CTIME'])
            else:
                ref = 1407000000
            reftime = int(ref) - 1407000000
        else:
            reftime = hdr[str(self.exptime_name)]
        return reftime

    def std_keywords(self, hdr):
        del hdr[('FILTER')]
        return

    #(1) It  modifies datasec and biassec of the header
    def modif_chan_header(self, data_dir, refname, chan, **kwargs):
        infile = os.path.join(data_dir, refname)
        hdulist = pyfits.open(infile)
        hdr = hdulist[0].header
        self.reftime = self.get_reftime((pyfits.open(infile))[0].header)
        if(self.device == 'Subaru'):
            self.correct_Sec_Subaru(hdr)
        else:
            self.correct_dataAndBiasSec(hdr, chan)
        image = hdulist[0].data
        hdu = pyfits.PrimaryHDU(data=image, header=hdr)
        hdu.writeto("%s"%infile, clobber=True)
        return infile

        #(1) it returns 2 frames : 0=flat1, 1=mask1,
        #(2) It also load mean and var of the overscan

        #(1) It copies a segment from a data repository to the current repository
        #(2) With the images from our testbench, ImageCopy doesn't work, I swith to imcopy
        #(3) It also modified datasec and biassec of the header
        #(4) Its gets the date of exposure in sec since 1.1.1970 minus 1407000000 sec.
    def local_copy_and_modif_chan_header(self, data_dir, refname, chan, method="ImageCopy", **kwargs):
        temporary_image = "temp.fits"
        os.system("rm -f %s" % temporary_image)
        infile = os.path.join(data_dir, refname)
        if(self.device == 'Subaru'):
            extension = "%s"%(infile)
        else:
            extension = "%s[%i]"%(infile, chan)

        if(method == "ImageCopy"):
            bf.ImageCopy(extension, temporary_image)
        if(method == "imcopy"):
            os.system("imcopy %s %s" % (str(extension), str(temporary_image)))

        hdulist = pyfits.open(temporary_image)
        hdr = hdulist[0].header

        self.reftime = self.get_reftime((pyfits.open(infile))[0].header)
        if(self.device == 'LapinDeJade'):
            print("Device is LapinDeJade")
            self.std_keywords(hdr)
            append_header_fromfile(hdr)
        if(self.device == 'Subaru'):
            self.correct_Sec_Subaru(hdr)
        else:
            self.correct_dataAndBiasSec(hdr, chan)
        image = hdulist[0].data
        hdu = pyfits.PrimaryHDU(data=image, header=hdr)

        hdu.writeto("%s"%temporary_image, clobber=True)
        return temporary_image

        #(1) it returns 2 frames : 0=flat1, 1=mask1,
        #(2) It also load mean and var of the overscan
    def img_and_masks(self, temporary_image, refname, chan):
        img_and_masks = []
        self.overscan_mean_and_var = []
        maskname = "mask_" + str(chan) + "_" + refname
        image = self.fitsfile2image(temporary_image)
        fitsimage = self.fitsfile2fitsimage(temporary_image)

        val = self.getsec(self.BIASSEC)
        self.overscan_mean_and_var.append(self.MeanVarValue(image, val)) #mean, var
        if(self.device == 'Subaru'):
            img_and_masks.append(self.mask_tearingSC2(fitsimage, name=maskname, amp=chan))
        elif(self.device == 'ITL'): # pas necesssaire
            #img_and_masks.append(self.mask_tearing(fitsimage, "image", name =  maskname)) #ok
            img_and_masks.append(self.mask_tearing(fitsimage, "fitsimage", name=maskname))
        else:
            img_and_masks.append(self.mask_tearing(fitsimage, "fitsimage", name=maskname))

        os.system("rm -f %s" % temporary_image)
        return img_and_masks

    #(1) It returns the mean and variance of the bias section
    def MeanVarValue(self, image, biassec):
        ccd = bf.ccdstat()
        bb = ccd.get_mean_var(image,
                              int(biassec[0]),
                              int(biassec[2]),
                              int(biassec[1]),
                              int(biassec[3]))
        return ccd.mean, ccd.var

    def fitsfile2image(self, fitsfile):
        return pc.Image(pc.FitsImage(fitsfile))

    def fitsfile2fitsimage(self, fitsfile):
        return pc.FitsImage(fitsfile)

    #(1) segment is either a fitsimage or an image
    #(2) It masks bad pixels using imageback.h algorithms
    def mask_tearing(self, fitsimage, TYPE, name="mask.fits", **kwargs):
        if(TYPE == ("image")):
            frame = fitsimage
        if(TYPE == ("fitsimage")):
            bb = bf.TotalIlluRegion(fitsimage)
            frame = fitsimage.Subimage(bb)
        mask = pc.FitsImage(str(name), frame.Nx(), frame.Ny())
        bf.mask_pixels(frame, mask, self.nsig)
        bf.PyCintex.gbl.convolve_mask(mask, self.bordersize, self.minpixcount)
        return frame, mask

    def mask_tearingSC2(self, fitsimage, name="mask.fits", amp=0, **kwargs):
        bb = bf.IlluRegion(fitsimage, amp)
        frame = fitsimage.Subimage(bb)
        mask = pc.FitsImage(str(name), frame.Nx(), frame.Ny())
        bf.mask_pixels(frame, mask, self.nsig)
        bf.PyCintex.gbl.convolve_mask(mask, self.bordersize, self.minpixcount)
        return frame, mask
