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
import toolbox as tb


if __name__ == "__main__":
    show_stuff = True
    files = sys.argv[1:]
    for file in files:

        dict, values, names = tb.readcat(file)
        columns = ('NUMBER',
                   'FLUX_ISO',
                   'FLUXERR_ISO',
                   'FLUX_ISOCOR',
                   'FLUXERR_ISOCOR',
                   'FLUX_APER1',
                   'FLUX_APER2',
                   'FLUX_APER3',
                   'FLUX_APER4',
                   'FLUX_APER5',
                   'FLUXERR_APER1',
                   'FLUXERR_APER2',
                   'FLUXERR_APER3',
                   'FLUXERR_APER4',
                   'FLUXERR_APER5',
                   'MAG_APER',
                   'MAGERR_APER',
                   'FLUX_AUTO',
                   'FLUXERR_AUTO',
                   'FLUX_BEST',
                   'FLUXERR_BEST',
                   'THRESHOLD',
                   'FLUX_MAX',
                   'X_IMAGE',
                   'Y_IMAGE',
                   'X2_IMAGE',
                   'Y2_IMAGE',
                   'XY_IMAGE',
                   'CXX_IMAGE',
                   'CYY_IMAGE',
                   'CXY_IMAGE',
                   'A_IMAGE',
                   'B_IMAGE',
                   'FLAGS')

        data = np.array(values)
        outcat = file
        print 'writing ', outcat
        print 'output entries : ', columns
        tb.DumpTuple(columns, zip(*data), outcat)
