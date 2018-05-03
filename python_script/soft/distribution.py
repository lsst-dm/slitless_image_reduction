#!/usr/bin/env python 
'''
return r realization of the parameter param foolowing the observed distribution 

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
'''

import os, sys, re
import numpy as np
import pylab as pl
import argparse


def readcat(cat):
    objs = [];   columns = []
    fp = open( cat, "r")
    lines = fp.readlines()
    for line in lines :
        if len(line.strip()) != 0 :
            if (line[0]=='#'):
                if (line[0:4] != "#end") :
                    column = re.sub('#|:|\\n','', line)
                    columns.append(column)
                continue
            else :
                objs.append(line.split())     
    fp.close()
    info  = np.rec.fromrecords(np.array(objs, dtype=float), names = columns)
    return info, columns



def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "return r realization of the parameter param foolowing the observed distribution "
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="return values from a given distribution")
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
 
    parser.add_argument('-n',"--number", type=int, 
	                help = "number of realization", 
	                default=1)   
    parser.add_argument('-a',"--param", type=str, 
	                help = "name of observable", 
	                default='O3')
  
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args    = grabargs()
    plot    = args.plot
    param   = args.param
    number  = args.number

    if param=='O3':     # 10 years dataset
        prefix = 'TO3'
    if param=='PWV':    # 10 years dataset
        prefix = 'TQV'
    if param=='AOD':    # this is only from 20 days in Oct. 2017
        prefix = 'AOD'
    if param=='CLOUDS': # 10 years dataset
        prefix = 'TAUTOT'

        
    print 'reading ',prefix+'_distrib.list'
    values, names = readcat(prefix+'_distrib.list')
    x   = values.field('bin')
    y   = values.field('value')
    val = np.random.choice(x, number, p=y)

    print param , val


    if plot :
        pl.bar(x,y, color='blue')
        pl.xlabel(param)
        n, bins, patches = pl.hist(val, len(x),
                                   normed=True,
                                   facecolor='g', alpha=0.75)
        pl.show()
