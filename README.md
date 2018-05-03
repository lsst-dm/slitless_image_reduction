# slitless_image_reduction
These scripts are used to reduced slitless images from CTIO 0.9-m telescope





Command line task exemple :

slitless_image_reduction/python_script/soft/process.py --aperture 40 --method 0 --order m+1 --disp --dark --cosmic --suffix aper40 -i data/f284.fits

Which uses classes from :
slitless_image_reduction/python_script/lsst/utils/