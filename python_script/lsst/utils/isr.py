#!/usr/bin/env python

from __future__ import absolute_import
#!/usr/bin/env python

from builtins import str
from builtins import range
from builtins import object
from . import telinst as instru
import os
import sys
import astropy.io.fits as pf
import numpy as np
from croaks import NTuple
from . import toolbox as tb
from . import reduceimage as ri
from . import extraction as ex
from . import flatfield as fd
from . import dispersion as disp
import logging


# ===========================================================
class isr(Exception):
    pass
# ===========================================================


''' A class to remove instrumental signatures'''


class isr(object):
    def __init__(self, spectrum, **kwargs):
        self.spectrum = spectrum

    ''' This is used to model the periodic pattern\
    imprint by the electronic on the background'''

    def ExtractProfile(self, footprint, mask, direction):
        Y = []
        if(direction == "y"):
            footprint = footprint.T
            mask = mask.T
        rows, cols = len(footprint), len(footprint[0])
        for i in range(0, rows):
            keep = footprint[i, :]
            remo = mask[i, :]
            keep = keep[remo == 0]
            if len(keep) is 0:
                flux = Y[-1]
            else:
                flux = np.median(keep)
            Y.append(flux)

        rows = np.arange(0, rows)
        file = str("background"+self.spectrum.order+".list")
        tb.DumpTuple(('rows', 'level'), (rows, Y), os.path.join(self.spectrum.out_dir, file))
        return Y

    '''Catch surrounding frame of footprint and masks from mask'''
    ''' Also return footprint's mask                           '''

    def Mask(self, inst, direction):
        which_amp = inst.IsInAmp(self.spectrum.y_start, self.spectrum.x_start,
                                 Trimmed=True)
        amp_frame = inst.TelAmpLimits(which_amp, Trimmed=True)
        masking = instru.telinst(self.spectrum.mask, verbose='')
        m = masking.Image(self.spectrum.mask)
        '''A surrounding frame to extract the periodic pattern in the background'''
        y_mask_start = max(self.spectrum.y_start-100, amp_frame[0]) # x<->y
        y_mask_end = min(self.spectrum.y_end+100, amp_frame[1]) # x<->y
        m = tb.Flip(m[y_mask_start:y_mask_end,
                      self.spectrum.x_start:self.spectrum.x_end],
                    flip=self.spectrum.flip)
        f = tb.Flip((inst.Image(self.spectrum.image))[y_mask_start:y_mask_end,
                                                      self.spectrum.x_start:self.spectrum.x_end],
                    flip=self.spectrum.flip)

        return self.ExtractProfile(f, m, direction)
