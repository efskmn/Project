#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #####################################
#
#
# Convert RA DEC coordinates given in degrees to galactic coordinates and plot. 
#
#
############################# BLOCK IMPORTS AND STUFF #################################

import numpy as np
from astropy.table import Table 
from astropy.coordinates import SkyCoord 
from astropy import units as u

##########################################################
# Copy and paste the ra(s) and dec(s) in degrees from the literature and 
# convert them to galactic coordinates to see which of them are in the field
##########################################################


myt = Table([ras,decs,names]) 
cc = SkyCoord(myt["col0"],  myt["col1"], frame='icrs', unit=(u.deg, u.deg)) 


fig, ax = pl.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(13,9), dpi=80, facecolor='w', edgecolo
r='k')                                                                                                             
for ind, i in enumerate(myt["col0"]):
    ax.plot(cc.galactic.l.value[ind], cc.galactic.b.value[ind],'ko',ms=4.)
    ax.text(cc.galactic.l.value[ind], cc.galactic.b.value[ind], myt["col2"][ind])
    ax.text(cc.galactic.l.value[ind]+0.2, cc.galactic.b.value[ind]+0.2, "("+str(cc.galactic.l.value[ind])+","+str(c
c.galactic.b.value[ind])+")")

#ax.axis([343,312,-2.5,+0.3])
ax.set_xlabel("longitude (deg)")
ax.set_ylabel("latitude (deg)")    
