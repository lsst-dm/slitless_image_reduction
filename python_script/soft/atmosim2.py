#!/usr/bin/env python

################################################################
#
# Script to simulate air transparency with LibRadTran
# With a pure absorbing atmosphere
# Here we vary PWV
#################################################################
import os
import re
import math
import numpy as np
import sys,getopt
from subprocess import Popen,PIPE, STDOUT, call

# LibRadTran installation directory
libradtranpath = os.getenv('LIBRADTRANDIR')+'/'


ZXX           = 'z'    # XX index for airmass z :   XX=int(10*z)
WVXX          = 'wv'   # XX index for PWV       :   XX=int(pwv*10)
OZXX          = 'oz'   # XX index for OZ        :   XX=int(oz/10)
LSST_Altitude = 2.241  # in k meters 
OBS_Altitude  = str(LSST_Altitude)
TOPDIR        = '/Users/augustinguyonnet/harvard/atmo_simu'



############################################################################
def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(f):
        os.makedirs(f)
#########################################################################


def usage():
    print "*******************************************************************"
    print sys.argv[0],' -z <airmass> -w <pwv> -o <oz>'
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)
    print "*******************************************************************"
    
    


def model(model = 'baseline', **kwargs):
    if (model =="baseline"):
        params = "atmosphere_file  %s\n\
        data_files_path /Users/augustinguyonnet/harvard/soft/libRadtran-2.0.1/data\n\
        source solar %s\n\
        mol_abs_param SBDART\n\
        phi 0.\n\
        phi0 0.\n\
        rte_solver disort\n\
        wavelength 300 1100\n\
        output_quantity transmittance\n\
        sza %s\n\
        aerosol_default\n\
        altitude 2.241\n\
        mol_modify H2O %s MM\n\
        mol_modify O3  %s DU\n\
        output_user lambda edir\n\
        quiet " % (libradtranpath+'data/atmmod/afglus.dat',
                   libradtranpath+'data/solar_flux/kurudz_1.0nm.dat',
                   str(sza),
                   pwv_str,
                   oz_str)
        params += "\naerosol_angstrom %s 0.0 "    %(str(0.0))
        
    if (model =="best"):
        params = "atmosphere_file  %s\n\
data_files_path /Users/augustinguyonnet/harvard/soft/libRadtran-2.0.1/data\n\
source solar %s\n\
mol_abs_param SBDART\n\
wavelength 300 1100\n\
output_quantity transmittance\n\
altitude 2.241\n\
mol_modify H2O 4. MM\n\
aerosol_default\n\
aerosol_angstrom 1. 0.02\n\
mol_modify O3 270 DU\n\
aerosol_modify tau550 set 0.02\n\
sza %s\n\
output_user lambda edir \n\
quiet " % (libradtranpath+'data/atmmod/afglus.dat',
           libradtranpath+'data/solar_flux/kurudz_1.0nm.dat',
           str(sza))
       


    if (model =="none"):
        params = "atmosphere_file  %s\n\
data_files_path /Users/augustinguyonnet/harvard/soft/libRadtran-2.0.1/data\n\
source solar %s\n\
mol_abs_param SBDART\n\
wavelength 300 1100\n\
output_quantity transmittance\n\
altitude 2.241\n\
mol_modify H2O 4. MM\n\
aerosol_default\n\
aerosol_angstrom 0. 0.0\n\
mol_modify O3 270 DU\n\
sza %s\n\
output_user lambda edir \n\
quiet " % (libradtranpath+'data/atmmod/afglus.dat',
           libradtranpath+'data/solar_flux/kurudz_1.0nm.dat',
           str(sza))
       
    return params

        
    
def main(argv):
    airmass_str=""
    pwv_str=""
    oz_str=""
    try:
        opts, args = getopt.getopt(argv,"hz:w:o:",["am=","pwv=","oz="])
    except getopt.GetoptError:
        print 'test.py -z <airmass> -w <pwv> -o <oz>'
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-z", "--airmass"):
            airmass_str = arg
        elif opt in ("-w", "--pwv"):
            pwv_str = arg
        elif opt in ("-o", "--oz"):
            oz_str = arg   
         
  
    if airmass_str=="":
        usage()
        sys.exit()

    if pwv_str=="":
        usage()
        sys.exit()

    if oz_str=="":
        usage()
        sys.exit()
	

    return float(airmass_str),float(pwv_str),float(oz_str)	
 
#-----------------------------------------------------------------------------
if __name__ == "__main__":

    airmass_num,pwv_num,oz_num=main(sys.argv[1:])
    
    
    print '--------------------------------------------'
    print ' 2) airmass = ', airmass_num
    print ' 2) pwv = ', pwv_num
    print ' 2) oz = ', oz_num
    print '--------------------------------------------'    
   
 
    # manage input and output directories 
    TOPDIR2=TOPDIR
    ensure_dir(TOPDIR2)
    INPUTDIR=TOPDIR2+'/'+'in'
    ensure_dir(INPUTDIR)
    OUTPUTDIR=TOPDIR2+'/'+'out'
    ensure_dir(OUTPUTDIR)

    #water vapor   
    pwv_val=pwv_num
    pwv_str=str(pwv_val)
    wvfileindex=int(100*pwv_val)


    # airmass
    airmass=airmass_num
    amfileindex=int(airmass_num*1000)

    # Ozone    
    oz_val=oz_num
    oz_str=str(oz_num)
    ozfileindex=int(oz_num/10.)


    # Convert airmass into zenith angle 
    sza=math.acos(1./airmass)*180./math.pi

    params = model(model='best')
    print params



    BaseFilename='RT_LS_pp_tp_sa_rt_z'+str(amfileindex)\
                  +'_'+WVXX+str(wvfileindex)\
                  +'_'+OZXX+str(ozfileindex)                   

    inp=os.path.join(INPUTDIR,  BaseFilename+'.INP')
    with open(inp,'w') as f:
        f.write(params)

    
    
    out=os.path.join(OUTPUTDIR, BaseFilename+'.OUT')

    verbose = True
    if verbose:
        print("Running uvspec with input file: ", inp)
        print("Output to file                : ", out)
 
    cmd = libradtranpath+'bin/uvspec '+  ' < ' + inp  +  ' > ' + out
    
    if verbose:
        print("uvspec cmd: ", cmd)
    p   = Popen(cmd,shell=True,stdout=PIPE)
    p.wait()
  

