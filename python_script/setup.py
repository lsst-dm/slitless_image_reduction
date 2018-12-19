#! /usr/bin/env python
# -*- coding: iso-8859-1 -*-
# ======================================================================
#
# Decam and LSST utils setup installation script
#
# setup.py: installation script.
#
# ======================================================================
#
# This script should be run to install the Decam python software
# To build the package, do:
#
#     python setup.py build
#
# To install it properly, do:
#
#     python setup.py install
#
# for a system installation (you need to be superuser).
#
# or
#     python setup.py install --prefix=<mydirectory>
#
# for a personal installation in the directory <mydirectory>. In
# this case you should modify the PYTHONPATH environment variable
# accordingly.
#
# ======================================================================

from distutils.core import setup, Extension

setup(name='decam_lsst',
      version='0.0',
      description='Decam software',
      author=["Augustin Guyonnet",
              ],
      author_email="guyonnet@lpnhe.in2p3.fr",
      #scripts = ['get_ROnoise.py'
      #           ],
      packages=['lsst.utils'],#['decam.utils', 'lsst.utils']
      ext_modules=[],
      data_files=[]#('merope-reduction/data', ['data/merope-events.rdb'])]
      )


# ======================================================================
