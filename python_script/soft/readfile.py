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
from __future__ import print_function

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
    show_stuff = True
    files = sys.argv[1:]
    for file in files:
        dict, values = readlist(file)

        if show_stuff:
            #print 'values : ', dict.values()
            #print 'keys : ', dict.keys()
            #print values.field('wg')
            print('Dict : ', list(dict.items()))
            print('Array : ', values.field)
            pl.plot(values.field('wg'), values.field('tatmo'))
            pl.show()
