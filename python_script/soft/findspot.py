#!/usr/bin/env python
'''
give a raw fitsimage, position of object in it,
-> return a spectrum

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''

import os
import sys
import re
import optparse
import subprocess
import numpy as np
import pylab as pl
import astropy.io.fits as pf


''' A class to find the brighter object fits image'''


class FindSpot(object):
    def __init__(self, Image, **kwargs):
        self.img = Image

    def ParseSexCat(self, cat):
        objs = []
        keys = []
        fp = open(cat, "r")
        lines = fp.readlines()
        for line in lines:
            words = line.split()
            if words[0] == "#":
                keys.append(words[2])
            else:
                objs.append(words)
        fp.close()
        info = np.rec.fromrecords(np.array(objs, dtype=float), names=keys)
        print 'data in ', cat, keys
        return info

    '''Use output catalog of Sextractor to find the brighter rounder object'''

    def Position(self):
        default = os.path.join('/Users/lpnhe/harvard/atmosx_git/data', 'defaultPS.sex')
        cmd = "sex -c %s %s  -CATALOG_NAME=se.list "%(default, self.img)
        print cmd
        os.system(cmd)
        info = self.ParseSexCat('se.list')

        '''Select flag=0, bigger of quite round objects, max best flux'''
        info = info[info.field('FLAGS') == 0]
        info = info[(np.sqrt(info.field('X2_IMAGE')) / np.sqrt(info.field('Y2_IMAGE')) < 1.2)
                    & (np.sqrt(info.field('X2_IMAGE')) / np.sqrt(info.field('Y2_IMAGE')) > .83)]
        info = info[(np.sqrt(info.field('X2_IMAGE')) * np.sqrt(info.field('Y2_IMAGE')))
                    == max(np.sqrt(info.field('X2_IMAGE')) * np.sqrt(info.field('Y2_IMAGE')))]
        print 'Object bkgd and flux in pixel maximum :',\
            info.field('BACKGROUND')[0], info.field('FLUX_MAX')[0]
        return info.field('X_IMAGE')[0], info.field('Y_IMAGE')[0]


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "find spot using Sextractor"

    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--img", type=str,
                      help="fits image",
                      default='plots/c7675g0008o44/image.fits')
    opts, args = parser.parse_args()
    return opts


if __name__ == "__main__":
    narg = len(sys.argv)
    args = grabargs()
    image = args.img
    print 'Find spot on image : ', image
    spot = FindSpot(image)
    x_object, y_object = spot.Position()
    print 'x_object, y_object = ', x_object, y_object
