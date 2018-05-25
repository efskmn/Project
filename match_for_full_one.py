#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #####################################
#
#
# Match stars from individual detectors with a consecutive detector.
#
#
############################# BLOCK IMPORTS AND STUFF #################################

import numpy as np 
import gc
from astropy import units as u 
from astropy.table import Table
from astropy.coordinates import SkyCoord
from glob import glob

#######################################################################################

def match_two_rawh5(ddict, name1, name2):
    '''
    In order to obtain ddict dictionary, first run "fulltable.py"
    The CALIB_LB_011???JK.rawh5 type of files should be given as input
    The two given inputs must be consecutive; either horizontailly or vertically
    Orientation of the detectors:
    1   5    9    13
    2   6   10    14
    3   7   11    15
    4   8   12    16
    '''
    
    d011 = ddict[name1]["data"]
    d021 = ddict[name2]["data"] 

    c011 = SkyCoord(d011["L"]*u.deg, d011["B"]*u.deg, frame='icrs') 
    c021 = SkyCoord(d021["L"]*u.deg, d021["B"]*u.deg, frame='icrs')
    id_, d2d, d3d = c011.match_to_catalog_sky(c021)
    matched = d021[id_][d2d.arcsec<0.2]
    second_data = np.delete(d021, id_, 0)
    
    return d011, second_data


#######################################################################################


fig, ax = pl.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(9,10), dpi=80, facecolor='w',edgecolor='k')
d091, d101 = match_two_rawh5(ddict, "CALIB_LB_091d075JK.rawh5", "CALIB_LB_101d075JK.rawh5")
ax.plot(d091["L"],d091["B"],'k.',ms=0.6)
ax.plot(d101["L"], d101["B"],'r.',ms=0.6)
pl.show()











