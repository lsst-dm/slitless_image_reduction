#!/usr/bin/env python

from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
import os
import sys
import re

# ===========================================================


class LSST_utils(Exception):
    pass
# ===========================================================


def segments():
    ret = ['%02d' % amp for amp in range(10, 18)]
    ret2 = ['%02d' % amp for amp in range(0, 8)]
    ret.extend(ret2[::-1])
    return (ret)


def amp():
    num = ['%d' % amp for amp in range(1, 17)]
    return num


def seg2amp(seg):
    seg = str(seg)
    if (len(list(str(seg))) != 2):
        seg = '0'+seg
    segm = segments()
    ampl = amp()
    c = dict(list(zip(segm, ampl)))
    return c[str(seg)]


def amp2seg(ampli):
    ampli = str(ampli)
    segm = segments()
    ampl = amp()
    c = dict(list(zip(ampl, segm)))
    return c[str(ampli)]

### 2 lists with pair images


def get_pairs(data_dir, ref_string):
    import glob
    set1 = ref_string + "*1.fits"
    set2 = ref_string + "*2.fits"

    set1 = sorted(glob.glob(os.path.join(data_dir, set1)))
    set2 = sorted(glob.glob(os.path.join(data_dir, set2)))

    if (len(set1) != len(set2)):
        print("list of flat pairs have unequal length [aborted] !")
        sys.exit()
    for i, j in zip(set1, set2):
        if str(i).replace("1.fits", "") != str(j).replace("2.fits", ""):
            print("Pairs are not identical [aborted]!")
            sys.exit()
    return set1, set2
