#!/usr/bin/env python
from __future__ import print_function
from builtins import range
import os
import sys
import numpy as np
import pylab
import pyfits
import optparse
import seg_output as seg
import poloka.core as pc
import pixel_utils as pu

code_name = "<mask_tearing.py> "


def read_option():
    usage = "usage: [%prog] [options]\n"
    usage += "mask bad pixels on bias frames and return mean signal and variance"

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

    parser.add_option("-o", "--out",
                      dest="outname",
                      type="string",
                      help="outfile name",
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

    option, args = parser.parse_args()
    return option

### 2 lists with pair images


def get_bias(data_dir, ref_string):
    import glob
    ref = ref_string + "*.fits"
    ref = sorted(glob.glob(os.path.join(data_dir, ref)))
    return ref


if __name__ == '__main__':
    print("from ~/python_script/soft/")
    option = read_option()
    data_dir = option.data_dir
    ref_string = option.refname
    outfile = option.outname
    all_amps = option.all_amps
    name1 = get_bias(data_dir, ref_string)
    ### mask parameters
    nsig = 4
    bordersize = 3
    minpixcount = 3
    rxx = pu.Rxx_flat(nsig, bordersize, minpixcount, biassec='[530:540,10:1990]', datasec='[20:515,10:1990]')
    if all_amps:
        init = 1
        fin = 17
    else:
        init = (option.AMP)
        fin = (option.AMP) + 1
    ### here it iterates on amplifiers
    nt = []
    for refname1 in (name1):
        refname1 = refname1.split('/')[-1]
        for chan in range(init, fin):
            covXX = []
            #( 1 ) it calculates the mean bias level and variance of each bias exposure
            temporary_image = rxx.local_copy_and_modif_chan_header(data_dir, refname1, chan)
            img_and_masks1 = np.array(rxx.img_and_masks(temporary_image, refname1, chan)).flatten()
            p1 = np.array(rxx.overscan_mean_and_var).flatten()[0]
            sp1 = np.sqrt(np.array(rxx.overscan_mean_and_var).flatten()[1])
            mean_im1 = pc.PyCintex.gbl.get_mean_Mask(img_and_masks1[0], img_and_masks1[1])
            var = pc.PyCintex.gbl.covariances_MaskOutliers(img_and_masks1[0], img_and_masks1[1], 0, 0)
            reject_pix = pc.PyCintex.gbl.nb_rejected_pixel(img_and_masks1[1])
            d = [chan, reject_pix, mean_im1, p1, sp1, var]
            nt.append(tuple(d))
            os.system("rm -f mask*.fits")

    # dump the result
    out = open(outfile, 'w')
    out.write("""# ampl : 
# reject_pix :
# mean1 : 
# ped1 : 
# sped1 : 
# var : 
# end 
""")

    print("nb de colonnes = ", len(nt[0]))
    for l in nt:
        if None in l:
            continue
        out.write('%i %i %f %f %f %f\n' % l)
    out.close()

    ### dump pair difference ###
    #head = pc.FitsHeader(os.path.join(data_dir, refname1))
    #pc.FitsImage("pair.fits", head, pair_diff)
    ###                      ###
