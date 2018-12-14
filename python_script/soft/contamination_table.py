#!/usr/bin/env python
'''
return plots of spectra

Author: Augustin Guyonnet
guyonnet@lpnhe.in2p3.fr
'''
from __future__ import print_function

import os
import sys
import re
import numpy as np
import pylab as pl
import toolbox as tb
import croaks
import scipy.interpolate as interp


if __name__ == "__main__":

    outfile = sys.argv[1]
    files = sys.argv[2:]
    values = []
    fig = pl.figure()
    for file in files:
        print(file)
        name = os.path.basename(file)
        value, keys = tb.readlist(file, ('wavelength', 'contamination'))
        values.append(value)
        pl.plot(value[0], value[1], label=name)

        pl.ylim(0.001, .02)
        pl.legend()

    pl.xlabel('wavelength (nm)')
    pl.ylabel('Contamination (fraction)')
    pl.title('Ronchi 200 - 2nd order light contamination')
    fig.savefig("contamination.pdf")
    pl.show()

    ### I guess I will fit a constant value
    tb.DumpTuple(('wavelength', 'contamination'), (second[:, 0], interp/second[:, 1]), outfile)
