#!/usr/bin/env python
'''
Plot simulated atmospheric transparency

Author: Augustin Guyonnet
aguyonnet@harvard.fas.edu
'''

import sys

import numpy as np
import matplotlib.pyplot as pl
from astropy.table import Table
from astropy.io import fits as pf
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation


def usage():
    print "Usage: plot_from_fitstable the simulated atmospheric transmission"
    print "Return a plot of the transmission"
    print


#  column 0 : count number
#  column 1 : aerosol value
#  column 2 : pwv value
#  column 3 : ozone value
#  column 6 : data start


# animation function. This is called sequentially
def animate(i, data, line):
    wght = data[0]
    out = data[1]
    x = wght[6:]
    y = out[i, 6:]
    ax.clear()
    pl.text(500, 0.3,
            'Aerosol %s' % (str(out[i, 1])), color='red', fontsize=10)
    pl.text(500, 0.2,
            'pwv %s' % (str(out[i, 2])), color='red', fontsize=10)
    pl.text(500, 0.1,
            'Ozone %s' % (str(out[i, 3])), color='red', fontsize=10)
    #line.set_data(x, y)
    line = pl.plot(x, y, label='time '+str(i))
    pl.xlabel('wavelength (nm)')
    pl.ylabel('Atmospheric transmission')
    pl.legend()
    print 'anim ', i
    return line,


def TableFromFits(Fits):
    fits = pf.open(Fits)
    tbdata = fits[0].data                  # the table is in extension 1

    fig1 = pl.figure()
    pl.xlabel('wavelength (nm)')
    pl.ylabel('Atmospheric transmission')
    pl.legend()

    fits.info()
    data = fits[0].data

    out = []
    for i, item in enumerate(data[1:, :]):
        if item[50] > 0.0:
            out.append(item)
            if plot:
                pl.clf()
                print i, item[:10]
                pl.plot(data[0, 6:], item[6:], label='time '+str(i))
                pl.text(500, 0.3,
                        'Aerosol %s' % (str(item[1])), color='red', fontsize=10)
                pl.text(500, 0.2,
                        'pwv %s' % (str(item[2])), color='red', fontsize=10)
                pl.text(500, 0.1,
                        'Ozone %s' % (str(item[3])), color='red', fontsize=10)
                pl.legend()
                pl.pause(0.25)

    if plot:
        pl.show()
    print data[0, :]

    out = np.array(out)
    return out, data[0, :]


if __name__ == "__main__":
    plot = False
    Writer = animation.writers['ffmpeg']
    name = 'ref'
    fits = sys.argv[1]
    print fits
    out, wght = TableFromFits(fits)
    print 'out ', out[0, :6]

    #pl.plot(wght[6:], out[0,6:])
    #pl.show()

    frames = len(out[:, 0])
    print 'nb of frames ', frames

    fig, ax = pl.subplots()
    line = ax.plot([], [])

    scat = [wght, out]
    anim = FuncAnimation(fig, animate, fargs=(scat, line),
                         frames=frames, interval=4000,
                         repeat=False)

    #pl.show()
    outname = "datachallenge"
    print 'writing movie : ', outname
    writer = Writer(fps=2, metadata=dict(artist='Me'), bitrate=10800)
    anim.save(outname+'.mp4', writer=writer)

    print 'pf'
