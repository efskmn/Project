#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #####################################
#
# Draw and select two masked regions for MS and RC. Obtain their histograms with the CMD plot.
# Example: in a working ipython window
# run #/roipoly_select_regions.py -dat "x.raw" -j 3 -k 5 -ms 0.3 -saveornot True
 
############################# BLOCK IMPORTS AND STUFF #################################

import sys
import os
import os.path
import argparse as argp
import pylab as pl
from roipoly import roipoly, roi
import numpy as np
import gc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from string_bool import str2bool

#######################################################################################
    
def chunks(l, n):                 
    for i in range(0, len(l), n): 
        yield list(l[i:i+n])      
    
    
############################# BLOCK READ INPUT ########################################

parser = argp.ArgumentParser(description='Plot and make two ROIs!')

parser.add_argument('-dat',help='dat: the ALLSTAR/DAOPHOT output to work with.',required=True)

parser.add_argument('-j',help='j: The column number(s) for j.',required=True, type=int)

parser.add_argument('-k',help='k: The column number(s) for k.',required=True, type=int)

parser.add_argument('-skip',help='skip: How many lines to skip while reading "dat".',required=False,type=int,default=3)

parser.add_argument('-ms', help='ms: marker size for the plot', required=False, type=float, default=0.3)

parser.add_argument('-saveornot1', help='saveornot1: want to save the strip plot?', required=False, type=str2bool,nargs='?', default=False)

parser.add_argument('-saveornot2', help='saveornot2: want to save the strip with histograms plot?', required=False, type=str2bool,nargs='?', default=False)

parser.add_argument('-mssim', help='mssim: Do you want to apply the MS strip with Cumulative Distribution Function?', required=False, type=str2bool,nargs='?', default=False)

parser.add_argument('-spop', help='spop: to which synthetic population?', required=False, type=str, nargs='?', default=False)


args = parser.parse_args()

#######################################################################################

#os.chdir("/net/nas3/popes/efsokmen/MW/DATA/VVV/DISK/d038")
os.chdir("/run/media/efsokmen/Elements/SIMULATIONS/SYNTHETIC_D_EJK/")

# Read the data
dat = np.genfromtxt(args.dat, usecols=(args.j,args.k),skip_header=3, invalid_raise=False)
print (np.shape(dat))

m = (dat[:,0] < 30) & (dat[:,0] > -30) & (dat[:,1] < 30) & (dat[:,1] > -30)

#Mask out for the data that is too large
X = dat[m,0]-dat[m,1]
Y = dat[m,1]

# show the image
figure, ax = pl.subplots(num=None,nrows=1, ncols=1, figsize=(11, 10), dpi=80, facecolor='w', edgecolor='k')
ax = pl.gca()
ax.invert_yaxis()
pl.plot(X ,Y, 'k.', ms=args.ms)
pl.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.6)
pl.title('1st ROI, left click: line segment, right click: close region')
ax.set_ylim(11,20.)
ax.set_xlim(-0.35,2.7)
ax.invert_yaxis()

# let user draw first ROI
ROI1 = roipoly(roicolor='r') #let user draw first ROI

# show the image with the first ROI
ax = pl.gca()
pl.plot(X ,Y, 'k.', ms=args.ms)
pl.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.6)
ax.set_ylim(11,20.)
ax.set_xlim(-0.35,2.7)
ax.invert_yaxis()
ROI1.displayROI()
pl.title('2nd ROI, left click: line segment, right click: close region')

# let user draw second ROI
ROI2 = roipoly(roicolor='b') #let user draw ROI

# show ROI masks
figure, ax = pl.subplots(num=None,nrows=1, ncols=1, figsize=(11, 10), dpi=80, facecolor='w', edgecolor='k')
ax.plot(X,Y,'k.', ms=args.ms)
ROI1.displayROI()
ROI2.displayROI()
pl.xlabel("J-K",fontsize=14)
pl.ylabel("K",fontsize=14)
pl.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.6)
ax.set_ylim(11,20.)
ax.set_xlim(-0.35,2.7)
ax = pl.gca()
ax.invert_yaxis()
pl.title('RC (red) and MS (blue) strips for {}'.format(args.dat),fontsize=14)
pl.show()

if args.saveornot1 is True:
    figure.savefig("roimasks"+args.dat[:-6]+".png",format='png')
    print ("Saved", "roimasks"+args.dat[:-6]+".png")
else:
    pass  


mask_roi1 = roi(ROI1, X , Y)
mask_roi2 = roi(ROI2, X, Y)

figure, ax = pl.subplots(num=None,nrows=1, ncols=1, figsize=(11, 10), dpi=80, facecolor='w', edgecolor='k')
ax.plot(X,Y,'k.', ms=args.ms)
ROI1.displayROI()
ROI2.displayROI()
pl.xlabel("J-K",fontsize=14)
pl.ylabel("K",fontsize=14)
pl.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)
ax = pl.gca()
ax.invert_yaxis()
ax.axis([0.0,2.7,19.8,12])
figure.text(0.20,0.95,'RC (red) and MS (blue) strips for {}'.format(args.dat),fontsize=14)

divider = make_axes_locatable(ax)
axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)


# now determine nice limits by hand:
binwidth = 0.05
xymax = max(np.max(np.abs(X)), np.max(np.abs(Y)))
lim = (int(xymax/binwidth) + 1)*binwidth

bins = np.arange(-lim, lim + binwidth, binwidth)
# Color histogram for RC
n_rc_1,bi_rc_1,patches_rc_1 = axHistx.hist(X[np.where(mask_roi1)], bins=bins, histtype="step", lw=1.2, color='red')
# Ks histogram for RC
n_rc_2,bi_rc_2,patches_rc_2 = axHisty.hist(Y[np.where(mask_roi1)], bins=bins, histtype="step", lw=1.2, color='red',orientation='horizontal')

# Color histogram for MS
n_ms_1,bi_ms_1,patches_ms_1 = axHistx.hist(X[np.where(mask_roi2)], bins=bins, histtype="step", lw=1.2, color='blue')
# Ks histogram for MS
n_ms_2,bi_ms_2,patches_ms_2 = axHisty.hist(Y[np.where(mask_roi2)], bins=bins, histtype="step", lw=1.2, color='blue',orientation='horizontal')

#Write the ratio of MS to RC in the histogram plots
axHistx.legend(['$MS/RC={}$'.format(round((np.max(n_ms_1)/np.max(n_rc_1)),1))], loc='best')
axHisty.legend(['$MS/RC={}$'.format(round((np.max(n_ms_2)/np.max(n_rc_2)),1))], loc='best')


# make some labels invisible
axHistx.xaxis.set_tick_params(labelbottom=False)
axHisty.yaxis.set_tick_params(labelleft=False)

# set the ticks of the y axis of vertical histogram and x ticks of horizontal histogram
axHistx.set_yticks(np.arange(0, np.max(n_ms_1), 400))
axHisty.set_xticks(np.arange(0, np.max(n_ms_2), 400))

# set grids
axHistx.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)
axHisty.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)

pl.show()

# save or not the final plot
if args.saveornot2 is True:
    figure.savefig("roimasks_hist_"+args.dat[:-6]+".png",format='png')
    print ("Saved", "roimasks_hist_"+args.dat[:-6]+".png")
else:
    pass 

color_hist_rc = [n_rc_1,bi_rc_1]
color_hist_ms = [n_ms_1,bi_ms_1]

k_hist_rc = [n_rc_2,bi_rc_2]
k_hist_ms = [n_ms_2,bi_ms_2]

# For the MS
sorted_jks = sorted(X[m][np.where(mask_roi2)])
minimum, maxsimum = np.min(sorted_jks), np.max(sorted_jks)


nbinsize=0.04
nbins = np.arange(minimum-0.5, maxsimum+0.5, nbinsize)

# Plot the MS histogram and shifted CDF
figur, axs = pl.subplots(nrows=1, ncols=2, sharex=False, sharey=False, figsize=(10,5), dpi=80, facecolor='w', edgecolor='k')
n,bi,pa = axs[0].hist(sorted_jks, bins=nbins, histtype="step", lw=1.2, color='black')

# Find how much the histogram has to be shifted in X to move the maximum of it to 0
binxy = np.vstack([bi[:-1]+nbinsize/2., n])
bin_m = (binxy[1] == np.max(binxy[1]))
shift_x = binxy[0][bin_m][0]

# Shifted histogram
nbins2 = nbins - shift_x

axs[0].hist(sorted_jks-shift_x*np.ones(len(sorted_jks)),bins=nbins2, histtype="step", lw=1.2, color='red', label="Shifted by {}".format(shift_x))
axs[1].step(np.array(sorted_jks), np.arange(X[m][np.where(mask_roi2)].size)/float(X[m][np.where(mask_roi2)].size), 'k--',lw=1.2)
axs[1].step(np.array(sorted_jks)-shift_x, np.arange(X[m][np.where(mask_roi2)].size)/float(X[m][np.where(mask_roi2)].size), 'r--',lw=1.2)
figur.text(0.5, 0.04, 'J-Ks', ha='center', fontsize=16)
axs[0].set_ylabel("Counts",fontsize=16)
axs[1].set_ylabel("CDF",fontsize=16)
figur.text(0.3,0.95, "MS strip of {}".format(args.dat[:-6]),fontsize=14)
axs[0].legend("best")
pl.show()


# Save the final plot and the corresponding .dat file or not. .dat file to be used for reddening and distance dispersions
if args.saveornot2 is True:
    figur.savefig("MS_strip_CDF_"+args.dat[:-6]+".png",format='png')
    print ("Saved", "MS_strip_CDF_"+args.dat[:-6]+".png")
    filoms = open("MS_strip_CDF_"+args.dat[:-6]+".dat", "w+")
    for i in np.arange(0,len(sorted_jks)):
        filoms.write( "%9.4f%9.4f\n" % ((np.array(sorted_jks)-shift_x)[i], (np.arange(X[m][np.where(mask_roi2)].size)/float(X[m][np.where(mask_roi2)].size))[i] ) )
    filoms.close()
    print ("Saved file", "MS_strip_CDF_"+args.dat[:-6]+".dat")
else:
    pass 
  

  
