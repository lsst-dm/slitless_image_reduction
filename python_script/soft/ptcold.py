#!/usr/bin/env python
'''
analyze ptc

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
from __future__ import print_function


import os
import sys
import astropy.io.fits as pf
import glob

#!/usr/bin/env python
import os
import sys
import numpy as np
import pylab
import pyfits
import optparse
import seg_output as seg
import poloka.bfstuff as bf
import poloka.core as pc
import pixel_utils as pu

code_name = "<mask_tearing.py> "


def read_option():
    usage = "usage: [%prog] [options]\n"
    usage += "mask bad pixels and return pixel correlation on a different of 2 flatfields with same illumination time"

    parser = optparse.OptionParser(usage=usage)

    parser.add_option("-d", "--dir",
                       dest="data_dir",
                       type="string",
                       help="repository origin of the images",
                       default="./")

    parser.add_option("-i", "--img",
                       dest="refname",
                       type="string",
                       help="image name",
                       default="./")

    parser.add_option("-v", "--verbose",
                       dest="verbose",
                       action="store_true",
                       help="decides to show or not commentaries on specific points of the code",
                       default=False)

    parser.add_option("--all_amps",
                       action="store_true",
                       dest="all_amps",
                       help="decides to run the script on all amps",
                       default=False)

    parser.add_option("-a", "--amp",
                       dest="AMP",
                       type="int",
                       help="run the script on amp",
                       default=None)

    parser.add_option("--renorm",
                       action="store_true",
                       dest="renorm",
                       help="It renormalizes img2 on img1",
                       default=False)

    option, args = parser.parse_args()

    return option

### 2 lists with pair images


def get_pairs(data_dir, ref_string):
    import glob
    set1 = ref_string + "*flat1*.fits"
    set2 = ref_string + "*flat2*.fits"
    set1 = sorted(glob.glob(os.path.join(data_dir, set1)))
    set2 = sorted(glob.glob(os.path.join(data_dir, set2)))

    if (len(set1) != len(set2)):
        print("list of flat pairs have inequal length !")
        sys.exit()
    for i, j in zip(set1, set2):
        print(i, j)
#        if str(i).replace("1.fits", "") != str(j).replace("2.fits", ""):
#            print "Pairs are not identical !"
#            sys.exit()
    return set1, set2


if __name__ == '__main__':
    option = read_option()
    data_dir = option.data_dir
    ref_string = option.refname
    all_amps = option.all_amps
    renorm = option.renorm
    name1, name2 = get_pairs(data_dir, ref_string)
    ### mask parameters
    nsig = 3;
    bordersize = 3;
    minpixcount = 3;
    rxx = pu.Rxx_flat(nsig, bordersize, minpixcount, biassec='[530:540,10:1990]', datasec='[38:510,20:1970]')
    if all_amps:
        init = 1
        fin = 17
    else:
        init = option.AMP
        fin = option.AMP + 1
    ### here it iterates on amplifiers
    nt = []
    for refname1, refname2 in zip(name1, name2):
        refname1 = refname1.split('/')[-1]
        refname2 = refname2.split('/')[-1]
        for chan in range(init, fin):
            covXX = []
            #( 1 ) it returns 4 frames : 0=flat1, 1=mask1, 2=flat2, 3=mask2
            #( 2 ) it calculates the mean flux of each flat, get the date of exposure in sec
                   ###  (overscan mean and var are calculated in (img_and_masks)
            #( 3 ) it makes the difference between the pairs
            #( 4 ) it masks the difference
            #( 5 ) it combines the 3 masks
        #( 6 ) it counts how many pixels are rejected on the segment
            #( 7 ) it returns the mean1, ped1, meam2, ped2, and covariances up to 4 pixels distance
            # IMAGE1
            temporary_image = rxx.local_copy_and_modif_chan_header(data_dir, refname1, chan)
            img_and_masks1 = np.array(rxx.img_and_masks(temporary_image, refname1, chan)).flatten()

            p1 = np.array(rxx.overscan_mean_and_var).flatten()[0]
            sp1 = np.sqrt(np.array(rxx.overscan_mean_and_var).flatten()[1])
            reftime1 = rxx.reftime
            # IMAGE2
            temporary_image = rxx.local_copy_and_modif_chan_header(data_dir, refname2, chan)
            img_and_masks2 = np.array(rxx.img_and_masks(temporary_image, refname2, chan)).flatten()
            p2 = np.array(rxx.overscan_mean_and_var).flatten()[0]
            sp2 = np.sqrt(np.array(rxx.overscan_mean_and_var).flatten()[1])
            reftime2 = rxx.reftime
            # PAIR
            mean_im1 = bf.PyCintex.gbl.get_mean_Mask(img_and_masks1[0], img_and_masks1[1])
            mean_im2 = bf.PyCintex.gbl.get_mean_Mask(img_and_masks2[0], img_and_masks2[1])
            ### Mean img2 is renormalized on mean im1, and then use (mean_im1 - p1) as flux ref
            ### If renormalized, careful before measuring the gain, 
            ### the coef 2 is no longer valide but becomes (1+alpha^2)  !!!         
            if renorm:
                alpha = (mean_im1 - p1) / (mean_im2 - p2)
            else :
                alpha = 1.0
            pair_diff       = (img_and_masks1[0] -p1) - (img_and_masks2[0] - p2)/(1/alpha)
            diff, diff_mask = rxx.mask_tearing(pair_diff, "image", "mask.fits")
            MASK            = img_and_masks1[1] + img_and_masks2[1] + diff_mask
            residu          = bf.PyCintex.gbl.get_mean_Mask(pair_diff, MASK )
            reject_pix      = bf.PyCintex.gbl.nb_rejected_pixel(MASK)
            for i in range(0, 5, 1):
                for j in range(0, 5, 1):
                    covXX.append(bf.PyCintex.gbl.covariances_MaskOutliers(diff, MASK, i,-j))

            d   = [chan, float(residu), int(reftime1), int(reftime2), alpha, reject_pix, mean_im1, p1, 
 sp1, mean_im2, p2,sp2]
            d   = d + list(np.array(covXX).flatten())
            nt.append(tuple(d))
            os.system("rm -f mask*.fits")
    
    # dump the result 
    if renorm:
        outname = 'cov_renorm.ntuple'
    else :
        outname = 'cov2.ntuple'
    out = open(outname, 'w')
    out.write("""# ampl : 
# residu :
# reftime1 : 
# reftime2 :
# alpha : 
# reject_pix :
# mean1 : 
# ped1 : 
# sped1 : 
# mean2 : 
# ped2 : 
# sped2 : 
# var : 
# cov01 : 
# cov02 : 
# cov03 : 
# cov04 : 
# cov10 : 
# cov11 : 
# cov12 : 
# cov13 : 
# cov14 : 
# cov20 : 
# cov21 : 
# cov22 : 
# cov23 : 
# cov24 : 
# cov30 : 
# cov31 : 
# cov32 : 
# cov33 : 
# cov34 : 
# cov40 : 
# cov41 : 
# cov42 : 
# cov43 : 
# cov44 : 
# end 
""")
    
    print("nb de colonnes = " , len(nt[0]))
    for l in nt:
        if None in l:
            continue
        out.write('%i %f %i %i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' % l)
    out.close()

