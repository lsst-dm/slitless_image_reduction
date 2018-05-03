#!/usr/bin/env python 
import os,re 
import numpy as np
import croaks
import logging
import astropy.io.fits as pf
from croaks import NTuple
import scipy.ndimage.filters as filt
import pylab as pl
import astropy.time
import dateutil.parser
from astropy.coordinates import Angle
from astropy import units as u




def ParallacticAngle(latitude, ha, dec):
    latitude = latitude.split( )
    latitude = float(latitude[0])- float(latitude[1])/60. - float(latitude[2])/3600.
    latitude = Angle(latitude, u.deg)
    latitude = latitude.radian
    ha       = Angle(ha, unit='hourangle')
    ha       = ha.radian
    dec      = Angle(dec, unit=u.deg)
    dec      = dec.radian
    parallactic_angle =\
            np.arctan( np.sin(ha) / ( np.cos(dec)*np.tan(latitude) - np.sin(dec)*np.cos(ha) ) )
    parallactic_angle = parallactic_angle*180/np.pi
    return parallactic_angle



def BkgdStat(footprint, mask):
    mean = np.mean( footprint[mask==0] )
    std  = np.std( footprint[mask==0] )
    return mean, std
    

            
'''subtract a longitudinal profile to the footprint'''
def SubtractProfile(footprint, Profile):
    for i in range(len(Profile)):
        footprint[: , i] = footprint[: , i] - Profile[i]
    return footprint
         

def Flip(frame, flip=False, **kwargs):
    if flip==True:
        frame = frame[:,::-1]
    return frame


def clean(list1, list2):
    clean        = np.array([list1, list2]).transpose()
    clean        = (clean[~np.isnan(clean).any(1)]).transpose()
    return clean[0,:], clean[1,:]

    
def flatten(seq):
    if not seq:
        return []
    elif isinstance(seq[0],list):
        return (flatten(seq[0])+flatten(seq[1:]))
    else:
        return [seq[0]]+flatten(seq[1:])


def readlist(inputfile, dic=[]):
    data = croaks.NTuple.fromtxt(inputfile)
  
    value =[]
    for name in dic:
        value.append(np.array(data[:][name]))
    return value, data.keys


def readfile(inputfile, dic=[]):
    data = croaks.NTuple.fromfile(inputfile)
    value =[]
    for name in dic:
        value.append(np.array(data[:][name]))
    return value


def readtxt(inputfile, dic=[]):
    data = croaks.NTuple.fromtxt(inputfile)
    value =[]
    for name in dic:
       value.append(np.array(data[:][name]))
    return value

   


def outKeys(dict):
    kk = list(dict.keys())
    kv = list(dict.values())
    outdict = {}
    for i,j in zip(kk,kv):
        outdict[i]= " ".join(j[:])
    return outdict


def readtuple(cat):
    objs = []; 
    fp = open( cat, "r")
    lines = fp.readlines()
    for line in lines :
        if (line[0]!='#'):
            objs.append(line.split())     
    fp.close()
    return objs


def readcat(cat):
    objs = [];   columns = []
    dict = {}
    fp = open( cat, "r")
    lines = fp.readlines()
    for line in lines :
        if len(line.strip()) != 0 :
            if (line[0]=='#'):
                if (line[0:4] != "#end") :
                    column = re.sub('#|:|\\n','', line)
                    columns.append(column)
                continue
            if line[0] == "@" :
                words = line[1:].split()
                dict[words[0]] = words[1:]
                continue
            else :
                objs.append(line.split())     
    fp.close()
    info  = np.rec.fromrecords(np.array(objs, dtype=float), names = columns)
    return dict, info, columns

def ParseSexCat(cat):
    objs = []; keys =[]
    fp = open( cat, "r")
    lines = fp.readlines()
    for line in lines :
        words = line.split()
        if words[0] == "#" :
            keys.append(words[2])
        else :
            objs.append(words )
    fp.close()
    info = np.rec.fromrecords(np.array(objs, dtype=float), names = keys)
    logging.debug('data in '+ str(cat) + str(keys))
    return info


def DumpTuple(names, list, file):
    out  = open(file, 'w')
    for i in names :
        out.write('# '+i + ' :' +'\n')
    out.write('#end'+ '\n')
    list = zip(*list)
    for i in list:
        out.write(' '.join(map("{}".format, i))+'\n')
    out.close()
    return


 
def DumpFile(keys, names, list, file):
    list = zip(*list)
    info = np.rec.fromrecords([i for i in list], names = names)
    info = info.view(NTuple)

    print type(keys)
    
    for key, value in keys.iteritems(): # to remove '' from items
        print key, value, type(value)
        if(type(value)==list):
            info.keys[key] = ' '.join(value)
        else:
            info.keys[key] = value
    info.totxt(file)
    return



def JD(image):
    date = pf.open(image)[0].header.get('DATE')
    dt   = dateutil.parser.parse(date)
    return  (astropy.time.Time(dt)).jd



def Filters(image):
    filters = pf.open(image)[0].header.get('FILTERS')
    return  filters


def Run(jd):
    run = 0
    if (jd < 2457460):
        run = 1
    if ((jd > 2457622) and (jd < 2457627)):
        run = 2
    if ((jd > 2457714) and (jd < 2457719)):
        run = 3
    if ((jd > 2457770.) and (jd < 2457774.)):
        run = 5
    if ((jd > 2457898.) and (jd < 2457920.)):
        run = 6
    if (jd > 2458022.):
        run = 7
    return run


'''look for known name reference'''
def STDname(name):
    print 'looking for :', name
    if (name.find('HD205905')>=0):
        target ='HD205905'
    elif ((name.lower()).find('hd14943')>=0):
        target = 'HD14943'
    elif (name.find('HD200654')>=0):
        target = 'HD200654'
    elif (name.find('HD185975')>=0):
        target = 'HD185975'
    elif ((name.lower()).find('mucol')>=0):
        target = 'MUCOL'
    elif (name.find('HR7950')>=0):
          target = 'HR7950'
    elif (name.find('HR9087')>=0):
          target = 'HR9087'
    elif ((name.find('HD108344_disp')>=0) or (name.find('HD108344')>=0)):
          target = 'HD108344'
    elif ((name.lower()).find('lamlep')>=0):
        target = 'lamlep'
    elif ((name.lower()).find('ksi02cet')>=0):
        target = 'ksi02cet'
    else :
        target = None
    print 'target = ', target
    return target



"""Should put here if prod on other than lpnlp250"""
def ProdInstall(prodpath = None, **kwargs):
    print 'prodpath  ', prodpath
    ''' Load environment variables'''
    #if os.environ.get('PROD_PREFIX') is None :
    if prodpath is None :
        execfile(os.environ.get('HOME')+"/harvard/proddir.py")
    else:
        execfile(os.path.join(prodpath,"proddir.py"))




'''Determination of the seeing at wght, from the sigma of gaussian fit orthogonal'''
'''to dispersion direction (return poly2 fit parameters)                         '''
def SeeingAtWght(w, psf_sigma):
    clean          = np.array([w, psf_sigma]).transpose()
    clean          = clean[(clean[:,0]>400) & (clean[:,0]<900)]
    seeing_at_wght = np.polyfit(clean[:,0], clean[:,1], deg=2)
    return seeing_at_wght 



'''Compute a seeing from the clump of star found from second moments in se.list'''
''' Does not work on ctio .9m images : the clump are the cosmics !             '''
def Seeing(cat, plot):
    data = ParseSexCat(cat)
    data = data[data.field('FLAGS')==0]
    '''flag = 0, S/N>5 cut, seeing < 5 pixels, ellipticity <20%'''  
    data = data[(data.field('FLUX_BEST')/data.field('FLUXERR_BEST')) >5]
    data = data[(data.field('X2_IMAGE')<5) & (data.field('Y2_IMAGE')<5)]
    data = data[(np.sqrt(data.field('X2_IMAGE'))>.6) & (np.sqrt(data.field('Y2_IMAGE'))>.6)]
    data = data[(np.sqrt(data.field('X2_IMAGE')) / np.sqrt(data.field('Y2_IMAGE'))<1.2)\
                & (np.sqrt(data.field('X2_IMAGE')) / np.sqrt(data.field('Y2_IMAGE'))>.83)]
    gmxx = data.field('X2_IMAGE')
    gmyy = data.field('Y2_IMAGE')
    gmxy = data.field('XY_IMAGE')
    # Select in moments space where the density is the highest=> Works fine
    # Always  check that the value here is realistic

    seeing_dispersion = 0.05
    density = 0
    for i,j,k in zip(gmxx, gmyy, gmxy):
        dist = np.sqrt(np.power(data.field('X2_IMAGE') - i,2)\
                       + np.power(data.field('Y2_IMAGE') - j,2))
        cutlist = data[dist < seeing_dispersion]
        number  = len(cutlist)
        if number > density :
            density = number
            seeing  = np.median(np.power(cutlist.field('X2_IMAGE') * cutlist.field('Y2_IMAGE')\
                               - cutlist.field('XY_IMAGE') * cutlist.field('XY_IMAGE'), .25))
            sx = np.std(np.sqrt(cutlist.field('X2_IMAGE')))/density # sx dispersion of the clump
            sy = np.std(np.sqrt(cutlist.field('Y2_IMAGE')))/density # sy disperison of the clump
            logging.debug( 'density seeing, sx, sy = '+ str(density)+' '+str(seeing)+' '+str(sx)+' '+str(sy))
            #print 'number, x, y, sx, sy, flux'
            #for i,j,k,l,m,n in zip(cutlist.field('NUMBER'),cutlist.field('X_IMAGE'),cutlist.field('Y_IMAGE'), np.sqrt(cutlist.field('X2_IMAGE')), np.sqrt(cutlist.field('Y2_IMAGE')), cutlist.field('FLUX_MAX')):
            #    print i,j,k,l,m,n
        
            print
            
    if plot is True :
        pl.plot(np.sqrt(data.field('X2_IMAGE')),  np.sqrt(data.field('Y2_IMAGE')), 'r^')
        pl.show()
    logging.info('IQ (pix.) = '+ str(seeing)+' sx=  '+str(sx)+' sy = '+str(sy))
    return seeing, sx, sy




''' Loading for atmospheric transmission '''
def AtmoT():
    if os.environ.get('ATMO_PREFIX') is None :
        execfile("/Users/lpnhe/harvard/proddir.py")
    data   = np.recfromtxt(os.path.join(os.environ['ATMO_PREFIX'],\
                                        'pwv00.0_aero0.00_alpha1.3.atm'))

    #fig = pl.figure()
    #pl.plot(data[:,0], data[:,1], color='b')
    #pl.xlabel('wavelength (nm)')
    #pl.ylabel('atmospheric transmission')
    #pl.title('CTIO - fiducial atmosheric transmission')
    #pl.show()
    #pl.legend()
    #fig.savefig('atmosphere.pdf')
    return data[:,0], data[:,1]



''' Loading mirror reflextivity '''
def MirrorR():
    if os.environ.get('INSTRU_PREFIX') is None :
        execfile("/Users/lpnhe/harvard/proddir.py")
    mirror = np.loadtxt(os.path.join(os.environ['INSTRU_PREFIX'], 'miror.dat'))
    mirror[:,0] = mirror[:,0]/10. # convert from A to nm
    return mirror[:,0], mirror[:,1]



''' Loading imager QE '''
''' Careful if using manuf. QE -> need to multiply by mirror reflectiviy'''
''' Should write a method telThroughput'''
def telescope_T(telescope):
    if os.environ.get('INSTRU_PREFIX') is None :
        execfile("/Users/lpnhe/harvard/proddir.py")
  
    if (telescope == 'open'):
        file = os.path.join(os.environ['INSTRU_PREFIX'],'Oct2017QE.list')
          #file = os.path.join(os.environ['INSTRU_PREFIX'],'QEOctober2017.txt') # not normalized by CBP throughput
        #file = os.path.join(os.environ['INSTRU_PREFIX'],'qecurve.txt')#the usual one
        #file = os.path.join(os.environ['INSTRU_PREFIX'],'qecurve_uv.txt')

    if (telescope == 'r'):
        file = os.path.join(os.environ['INSTRU_PREFIX'],'r_band2.list')
        
    if (telescope == 'rg715'): 
        #file = os.path.join(os.environ['INSTRU_PREFIX'],'r400_obf.list')# bad measure
        file = os.path.join(os.environ['INSTRU_PREFIX'],'qecurve.txt')
        
    qe  = np.loadtxt(file)
    #fig = pl.figure()
    #norm   = max(qe[:,1])
    #signal = qe[:,1]/norm
    #pl.plot(qe[:,0], signal, color='b')
    #pl.xlabel('wavelength (nm)')
    #pl.ylabel('Telescope throughput (norm)')
    #pl.title('CTIO 0.9m - CBP calibration')
    #pl.show()
    #pl.legend()
    #fig.savefig('tel_throughput.pdf')
    print 'Telescope throughput from : ', file
    return   qe[:,0], qe[:,1]


def ccd_qe():
    file = os.path.join(os.environ['INSTRU_PREFIX'],'qecurve.txt')
    qe  = np.loadtxt(file)
    print 'QE from : ', file
    return   qe[:,0], qe[:,1]


'''Needed only if using CCD QE from manufacturer curve''' 
def QEMirror():
    wg, qe       = ccd_qe()
    wg2, R       = MirrorR()
    R            = interp.griddata(wg2, R, wg)
    qeR          = qe * R * R
    return wg, qeR



def BkgdStatClip(data, control_plot):
    lenX  = len(data)
    lenY  = len(data[0])
    c, low, up = stats.sigmaclip(data, 5, 2)
    print 'Lower and upper value included in bkgd : ', low, up
    mean  = np.mean(c)
    sigma = np.std(c)
    print 'Bkgd mean and sigma = ', mean, sigma
    if (control_plot is True):
        check = np.zeros([lenX, lenY])
        check = np.where(((data<up) & (data>low)), 1, 0)    
        pl.pcolor(check, cmap='hot')
        pl.colorbar()
        pl.show()
    return mean, sigma
    




def Closest(myList, myNumber):
    val = min(myList, key=lambda x:abs(x-myNumber))
    for i,j in enumerate(myList) :
        if (j==val):
            item = i
            break
    return item, val


''' run sectractor on image to extract a background map and a se.list'''
def RunSEx(image, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    home = os.environ['HOME']
    os.environ['SEX_PREFIX'] = os.path.join(home, 'harvard/soft/sextractor-2.5.0/config')
    default = os.path.join(os.environ['SEX_PREFIX'], 'default.sex')
    cmd = "sex -c %s %s  -CATALOG_NAME=%s"%(default, image,
                                            os.path.join(outdir , 'se.list'))
    print cmd
    os.system(cmd)
    os.system('mv segmentation.fits %s' %(outdir))
    return
     
