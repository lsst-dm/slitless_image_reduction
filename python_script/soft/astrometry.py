#!/usr/bin/env python
'''
just this :
http://docs.astropy.org/en/stable/wcs/

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''
#from __future__ import division, print_function

import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import argparse
import dateutil.parser

from astropy import time
from astropy import wcs
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import ICRS
from astropy.coordinates import FK5
from astropy.coordinates import match_coordinates_sky

'''
If matches are bad, maybe use
 sip_foc2pix()
to account for wcs distorsions
'''


test = [['TYC_6391-1496-1', 337.357346, -20.947772]]
references = [['cd-34241', '00  41  46.9212', '-33  39  08.430'],
              ['hip3142', '00  39  57.8235', '-33  57  41.218'],
              ['hr9087', '00  01  49.44853', '-03  01  39.0095'],
              ['cd-28595', '01  54  50.2715', '-27  28  35.747'],
              ['hip8417', '01  48  35.9288', '-26  15  13.147'],
              ['hd14943', '02  22  54.6751', '-51  05  31.660'],
              ['ksi02cet', '02  28  09.54266', '+08  27  36.2007'],
              ['cpd-69177', '03  10  31.02311', '-68  36  03.3866'],
              ['hip15968', '03  25  36.25177', '-69  20  11.1697'],
              ['hip17819', '03  48  46.8344', ' -40  23  57.216'],
              ['lp995-86', '03  48  22.61', '-39  08  36.9'],
              ['hip19796', '04  14  34.3409', '+10  42  04.992'],
              ['hr1544', '04  50  36.72298', '+08  54  00.6493'],
              ['hz2', '04  12  43.551', ' +11  51  48.75'],
              ['gd71', '05  52  27.614', '+15  53  13.75'],
              ['hip26382', '05  37  03.73543', '+17  02  25.1776'],
              ['lamlep', '05  19  34.52405', '-13  10  36.4408'],
              ['mucol', '05  45  59.89496', '-32  18  23.1630'],
              ['eg131', '19  20  34.9231', '-07  40  00.065'],
              ['eggr150', '21  52  25.38190', '+02  23  19.5405'],
              ['hip107120', '21  41  56.4652', '+00  20  45.778'],
              ['hip111449', '22  34  41.63641', '-20  42  29.5779'],
              ['hip112542', '22  47  42.76932', '-14  03  23.1409'],
              ['lp877-23', '22  52  41.032', '-20  35  32.89'],
              ['ngc7293', '22  29  38.541', '-20  50  13.64'],
              ['feige110', '23  19  58.39814', '-05  09  56.1604'],
              ['hip113896', '23  03  57.2734', '-04  47  41.502']]


def buildTargetCoords(references):
    reference = []
    for i in references:
        obj = i[0]
        ra = i[1]
        dec = i[2]
        ra, dec = sexa2deg(ra, dec)
        reference.append([obj, ra, dec])
    return reference


def sexa2deg(ra, dec):
    print 'was ra, dec : ', ra, dec
    ra = Angle(ra, unit='hourangle').degree
    dec = Angle(dec, unit=u.deg).degree
    print 'is in degree :', ra, dec

    return ra, dec


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "pixel to sky"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract spectrum from raw image")

    parser.add_argument('-v', "--output", type=str,
                        help="output catalog",
                        default=None)
    parser.add_argument('-i', "--reps", nargs='+', type=str,
                        help="list of reps to be processed",
                        default=None)
    parser.add_argument('-x', "--extern", type=str,
                        help="external wcs",
                        default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = grabargs()
    reps = args.reps
    extern = args.extern
    output = args.output
    if output:
        outcat = output
    else:
        outcat = 'catalog.list'
    columns = []
    data = []
    name = 0
    date = 0
    band = 0
    ra = 0
    dec = 0
    exptime = 0

    targets = buildTargetCoords(references)
    print targets

    '''Searching for matches'''
    for item, rep in enumerate(reps):
        file = os.path.join(rep, 'se.list')
        img = os.path.join(rep, 'calibrated.fits')
        print 'reading : ', file
        print 'opening img : ', img
        dict, objects, names = tb.readcat(file)

        if extern:
            hdulist = fits.open(extern)
        else:
            hdulist = fits.open(img)
            date = hdulist[0].header.get('DATE')
            date = dateutil.parser.parse(date)
            date = (time.Time(date)).jd
            band = hdulist[0].header.get('FILTERS').replace(" ", "_")
            exptime = hdulist[0].header.get('EXPTIME')
            airmass = hdulist[0].header.get('AIRMASS')
            outhum = hdulist[0].header.get('OUTHUM')
            outpress = hdulist[0].header.get('OUTPRESS')
            wndspeed = hdulist[0].header.get('WNDSPEED')

        '''Maybe match keyword name and reference name...unreliable'''
        head = (fits.open(img))[0].header
        obj = head.get('OBJECT')
        print 'OBJECT keyword : ', obj

        '''Convert sources coordinates from pixel to sky ra-dec'''
        w = wcs.WCS(hdulist[0].header)
        sources_ra = []
        sources_dec = []
        sources_nb = []
        enum = 0
        for i, j, k, l in zip(objects.field('X_IMAGE'),
                              objects.field('Y_IMAGE'),
                              objects.field('FLAGS'),
                              objects.field('NUMBER')):
            if(k == 0): # the one we look for is sometimes saturated
                pixcrd = np.array([[i, j]], np.float_)
                world = w.wcs_pix2world(pixcrd, 1)
                #print enum , 'Pix i,j :', i,j , ' is at ra, dec :', world[0][0], world[0][1]
                sources_ra.append(world[0][0])
                sources_dec.append(world[0][1])
                sources_nb.append(l)
                enum += 1

        catalog = FK5(Angle(sources_ra, u.degree), Angle(sources_dec, u.degree))

        for reference in targets:
            looking = reference[0]
            looking_ra = Angle(reference[1] * u.deg).degree
            looking_dec = Angle(reference[2] * u.deg).degree
            pixcrd = np.array([[reference[1], reference[2]]], np.float_)
            world = w.wcs_world2pix(pixcrd, 1)

            if((world[0][0] < 0) or (world[0][0] > 2000)
               or (world[0][1] < 0) or (world[0][1] > 2000)
               or (np.isnan(world[0][0]) == True) or (np.isnan(world[0][1]) == True)):
                continue
            print 'Continue only if is on focal plan'
            print 'Obj : ', looking, looking_ra, looking_dec
            print 'Is expected at pixels : ', world[0][0], world[0][1]

            '''look for the closest object '''
            c = FK5(Angle(looking_ra, u.degree), Angle(looking_dec, u.degree))
            idx, d2d, d3d = match_coordinates_sky(c, catalog)
            print idx, ' distance ~ ', (d2d.degree)[0]*3600 / 0.8,  #mode rebin 2*2
            print ' object is at number : ', sources_nb[idx]

            found = objects[objects.field('NUMBER') == sources_nb[idx]]
            print 'Found : ', found, ' => ', found.field('X_IMAGE'), found.field('Y_IMAGE')

            pixcrd = np.array([[found.field('X_IMAGE')[0], found.field('Y_IMAGE')[0]]], np.float_)
            world = w.wcs_pix2world(pixcrd, 1)
            ra = world[0][0]
            dec = world[0][1]
            name = reference[0]
            data.append(list([name, date, exptime, band, ra, dec, airmass,
                              outhum, outpress, wndspeed]) + list(found[0]))
            if item == 0:
                columns = ('name', 'date', 'exptime', 'band', 'ra', 'dec', 'airmass',
                           'outhum', 'outpress', 'wndspeed') + found.dtype.names
                print 'output entries : ', columns

    print 'writing ', outcat
    tb.DumpTuple(columns, zip(*data), outcat)
