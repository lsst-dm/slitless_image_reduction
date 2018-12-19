#!/usr/bin/env python
'''
replace keywords

Author: Augustin Guyonnet
guyonnet@fas.harvard.edu
'''

from builtins import zip
from builtins import str
import os
import sys
import astropy.io.fits as pf
import argparse


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "change keyword value"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="imgs, kw_names, kw_vals")
    parser.add_argument('-i', "--imgs", nargs='+', type=str,
                        help="fits image",
                        default=None)
    parser.add_argument('-k', "--keyword", nargs='+', type=str,
                        help="keyword names",
                        default=None)
    parser.add_argument('-v', "--value", nargs='+', type=str,
                        help="new values",
                        default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = grabargs()
    images = args.imgs
    keyword = args.keyword
    value = args.value

    for img in images:
        for k, v in zip(keyword, value):
            data, hdr = pf.getdata(img, header=True)
            old = hdr.get(k)
            comment = str('change from '+old)
            hdr.set(k, v, comment)
    pf.writeto(img, data, hdr, clobber=True)
