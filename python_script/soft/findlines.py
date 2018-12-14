#!/usr/bin/env python
'''
return plots of spectra

Author: Augustin Guyonnet
aguyonnet@fas.harvard.edu
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
import scipy.ndimage.filters as filt
import argparse


def Noise(data, SN):
    noise = []
    for i in data:
        noise.append(np.random.normal(i, (i/float(SN))))
    return noise


def Count(line):
    lines = []
    for i in line:
        for j in i:
            lines.append(j)
    lines = list(set(lines))
    return len(lines)


def FillSmooth2(ref, cleanedref):
    newref = interp.griddata(cleanedref[:, 0], cleanedref[:, 1], ref[:, 0])
    out = np.array([ref[:, 0], newref]).transpose()
    return out


def FillSmooth(ref, cleanedref, smooth):
    #bkgd   = interp.griddata(cleanedref[:,0], smooth, ref[:,0])
    lines = [i for i in ref[:, 0] if i not in cleanedref[:, 0]]
    fill = [(i, j) for i, j in zip(ref[:, 0], smooth) if i in lines]
    newref = np.concatenate((cleanedref, fill), axis=0)
    newref = newref[newref[:, 0].argsort()]
    return newref


def Closest(myList, myNumber):

    val = min(myList, key=lambda x: abs(x-myNumber))
    for i, j in enumerate(myList):
        if (j == val):
            item = i
            break
    return item, val


'''-------------------------------'''


class Lines(object):
    def __init__(self, data):
        self.data = data
        self.segment = []

    ''' Filter the spectrum using a median and iteratively separated       '''
    ''' the lines from the continuum based on the rms of spectrum-smoothing'''

    def Continuum(self):
        lendata = len(self.data[:, 0])
        size = lendata/20
        spec = self.data
        Iter = 1000
        line = []
        print('median over :', size)
        for i in range(Iter+1):
            smooth = filt.median_filter(spec[:, 1], size=size, mode='nearest')
            error = np.std(spec[:, 1]-smooth)
            if (error == 0.0):
                print('residuals of continuum reached 0. Leave iteration')
                data_continuum = spec
                break
            data_continuum = spec[(spec[:, 1]-smooth) >= (-3*error)]
            remove = list(spec[:, 0][(spec[:, 1]-smooth) < (-3*error)])
            line.append(remove)
            count = Count(line)
            lenremove = len(remove)
            size = int((lendata-len(line))/20)
            print('(iter, size,  error, removed) = ', i, size, error, count)

            if(lenremove == 0):
                print('No more point removed after iteration ', i, ' -> break')
                break
            '''replace the lines data point by the smooth coninuum estimation'''

            spec = FillSmooth2(spec, data_continuum)

            if (i == Iter):
                print('iter limit has been reached')
                print('return smoothed continuum out of last data selection')
                break

        lines = []
        for i in line:
            for j in i:
                lines.append(j)
        lines = list(set(lines))
        self.error = error
        return spec, lines

    def LinesLocation(self, continuum, lines):
        bkgd = [i for i in self.data[:, 0] if i not in lines]
        lines = np.array([lines, np.ones(len(lines))*float('nan')]).transpose()
        bkgd = np.array([bkgd, np.ones(len(bkgd))]).transpose()
        ref = np.concatenate((lines, bkgd), axis=0)
        ref = ref[ref[:, 0].argsort()]
        return ref[:, 1]

    def EquivalentWidth(self, bkgd, lines_loc):
        wg_lines = []
        flux_lines = []
        sflux_lines = []
        wg_sum = []
        flux_sum = []
        n = 0
        for i, j, k, l in zip(self.data[:, 0], bkgd, self.data[:, 1], lines_loc):
            if (np.isnan(l) and (j != k)):              # there is a line
                flux_sum.append(k-j)
                wg_sum.append(i)
                n += 1
            elif flux_sum:
                wg_lines.append(np.mean(wg_sum))
                flux_lines.append(sum(flux_sum))
                sflux_lines.append(np.sqrt(n) * self.error)
                flux_sum = []
                wg_sum = []
                n = 0
        print('lines position and flux : ')
        kept = []
        for i, j, k in zip(wg_lines, flux_lines, sflux_lines):
            if(np.abs(j/k) >= 2):
                kept.append(i)
                print(i, j, k)
        self.RefineLines(lines_loc, kept)
        return wg_lines, flux_lines, sflux_lines

    def RefineLines(self, bkgd, kept):
        wgth = []
        flux = []
        for i in kept:
            w, f = self.Flux(bkgd, i)
            wgth.append(w)
            flux.append(f)
        return wgth, flux

    def Flux(self, bkgd, line):
        dataC = np.array([self.data[:, 0], self.data[:, 1], bkgd]).transpose()
        lmin = max(dataC[:, 0][(dataC[:, 0] < line) & ((dataC[:, 2] == 1.))])
        min_point = self.Point(dataC, lmin)
        lmax = min(dataC[:, 0][(dataC[:, 0] > line) & ((dataC[:, 2] == 1.))])
        max_point = self.Point(dataC, lmax)
        position, flux, sflux = self.Sum(dataC, min_point, max_point)
        print("position, flux, sflux = ", position, flux, sflux)
        return position, flux

    def Point(self, dataC, point):
        # Instead of a mean, I should fit a spline and retrieve the highest derivative
        item, val = Closest(dataC[:, 0], point)
        selection = dataC[(item-2):(item+3), :]
        selection = selection[selection[:, 2] == 1.]
        p_pos = np.mean(selection[:, 0])
        p_val = np.mean(selection[:, 1])
        return p_pos, p_val

    def Sum(self, dataC, min_point, max_point):
        dataC = dataC[(dataC[:, 0] >= min_point[0]) & (dataC[:, 0] <= max_point[0])]
        linear_coef = np.polyfit([min_point[0], max_point[0]],
                                 [min_point[1], max_point[1]], deg=1)
        bkgd = np.polyval(linear_coef, dataC[:, 0])
        '''appending segment : dataC[:,0], bkgd'''
        self.segment.append([dataC[:, 0], bkgd])
        result = dataC[:, 1] - bkgd
        sflux = np.sqrt(len(result)) * self.error
        position = np.average(dataC[:, 0], weights=result)
        result = np.sum(result)
        return position, result, sflux

    def Plot(self, continuum):
        fig = pl.figure(0)
        pl.plot(self.data[:, 0], self.data[:, 1], color='black', label='Lines')
        pl.plot(continuum[:, 0], continuum[:, 1], color='crimson', label='Continuum')

        for i in self.segment:
            pl.plot(i[0], i[1], color='blue', linewidth=2)
        pl.title('Lines extraction')
        pl.xlabel('Wavelength (nm)')
        pl.ylabel('Spectrum Profile (arbitrary units)')
        pl.legend(loc='upper right')
        fig.savefig("lines2.pdf")
        return


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "extract lines from a spectrum"

    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract lines from a specturm")
    parser.add_argument('-f', "--file", type=str,
                        help="spectrum file",
                        default=None)
    parser.add_argument('-n', "--noise", type=int,
                        help="add Noise to spectrum with value from S/N",
                        default=None)
    parser.add_argument('-c', "--columns", nargs='+', type=int,
                        help="wavelength and flux columns in file",
                        default=[0, 1])
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = grabargs()
    if(args.columns):
        col0 = int(args.columns[0])
        col1 = int(args.columns[1])
    else:
        col0 = 0
        col1 = 1

    ''' input spectrum '''
    data = np.loadtxt(args.file, comments=["@", '#'], usecols=(col0, col1))

    ''' add noise'''
    if args.noise:
        SN = args.noise
        data[:, 1] = Noise(data[:, 1], SN)

    '''Declare class      '''
    ''' extract continuum '''
    li = Lines(data)
    continuum, lines = li.Continuum()

    ''' pinpoint line features in data'''
    lines_loc = li.LinesLocation(continuum, lines)

    ''' measure flux and position for each line '''
    wg_lines, flux_lines, sflux_lines =\
        li.EquivalentWidth(continuum[:, 1], lines_loc)

    ''' checks '''
    li.Plot(continuum)

    fig = pl.figure(1)
    pl.errorbar(wg_lines, flux_lines, yerr=sflux_lines,
                color='red', marker='^', linestyle='None')
    pl.xlabel('position (nm)')
    pl.title('Lines')
    #pl.legend()
    pl.show()
