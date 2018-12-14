#!/usr/bin/env python
'''
read and parse files with :
@ keys
# column_1
.
# column_n
#end

recarray doc :
https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.recarray.html

'''

import os
import sys
import re
import numpy as np
import pylab as pl


def readlist(cat):
    objs = []
    columns = []
    dict = {}
    fp = open(cat, "r")
    lines = fp.readlines()
    for line in lines:
        if len(line.strip()) != 0:
            if (line[0] == '#'):
                if (line[0:4] != "#end"):
                    column = re.sub('#|:|\\n', '', line)
                    columns.append(column)
                continue
            if line[0] == "@":
                words = line[1:].split()
                dict[words[0]] = words[1:]
                continue
            else:
                objs.append(line.split())
    fp.close()
    info = np.rec.fromrecords(np.array(objs, dtype=float), names=columns)
    return dict, info


if __name__ == "__main__":
    pwv = np.arange(0., 20., 0.5)
    aero = np.arange(0., 0.1, 0.02)
    ff = 1.3

    for p in pwv:
        for a in aero:
            params = "atmosphere_file ../data/atmmod/afglus.dat\n\
    source solar ../data/solar_flux/kurudz_1.0nm.dat\n\
            mol_abs_param SBDART\n\
            rte_solver disort\n\
            wavelength 300 1100\n\
            output_quantity transmittance\n\
            pressure 630\n\
            sza 62\n\
            aerosol_default\n\
            aerosol_angstrom %3.1f 0.02\n\
            aerosol_modify tau550 set %3.2f\n\
            altitude 2.241\n\
            mol_modify H2O %04.1f MM\n\
            mol_modify O3 270 DU\n\
            output_user lambda edir\n\
            quiet" % (ff, a, p)

            fname = 'pwv%04.1f_aero%3.2f_alpha%3.1f' % (p, a, ff)
            with open(fname, 'w') as f:
                f.write(params)
