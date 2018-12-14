#!/usr/bin/env python

import os
import sys
import re

# ===========================================================


class LSST_utils(Exception):
    pass
# ===========================================================


def segments():
    ret = ['%02d' % amp for amp in xrange(10, 18)]
    ret2 = ['%02d' % amp for amp in xrange(0, 8)]
    ret.extend(ret2[::-1])
    return (ret)


def amp():
    num = ['%d' % amp for amp in xrange(1, 17)]
    return num


def seg2amp(seg):
    seg = str(seg)
    if (len(list(str(seg))) != 2):
        seg = '0'+seg
    segm = segments()
    ampl = amp()
    c = dict(zip(segm, ampl))
    return c[str(seg)]


def amp2seg(ampli):
    ampli = str(ampli)
    segm = segments()
    ampl = amp()
    c = dict(zip(ampl, segm))
    return c[str(ampli)]

### 2 lists with pair images


def get_pairs(data_dir, ref_string):
    import glob
    set1 = ref_string + "*1.fits"
    set2 = ref_string + "*2.fits"

    set1 = sorted(glob.glob(os.path.join(data_dir, set1)))
    set2 = sorted(glob.glob(os.path.join(data_dir, set2)))

    if (len(set1) != len(set2)):
        print "list of flat pairs have unequal length [aborted] !"
        sys.exit()
    for i, j in zip(set1, set2):
        if str(i).replace("1.fits", "") != str(j).replace("2.fits", ""):
            print "Pairs are not identical [aborted]!"
            sys.exit()
    return set1, set2
