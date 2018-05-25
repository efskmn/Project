#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #############################
#
# Generic python script description.
# 
############################# BLOCK IMPORTS AND STUFF #########################
import sys
import os
import argparse as argp
import numpy as np
import os.path
import random as rdm    

### Pretty!
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
############################# BLOCK READ INPUT #############################

parser = argp.ArgumentParser(description='Resample, shuffle and divide!')

parser.add_argument('-xy',help='xy: DAOPHOT/ALLFRAME/ALLSTAR related file, with XY coords in 2nd and 3rd cols. ',required=True)

parser.add_argument('-sint',help='sint: File you need to add XY columns to.',required=True)

parser.add_argument('-skip',help='skip: How many lines to skip while reading "xy".',required=False,type=int,default=3)

parser.add_argument('-nsplit',help='nsplit: How many files you want to split "sint" into.',type=int,required=False,default=2)

args = parser.parse_args()

# This to ensure that the files actually exist...
if not os.path.isfile(args.sint):
  sys.exit(bcolors.FAIL + 'File ' + args.sint + ' not found!' + bcolors.ENDC)

if not os.path.isfile(args.xy):
  sys.exit(bcolors.FAIL + 'File ' + args.xy + ' not found!' + bcolors.ENDC)

# Read the XY from a RAW file.
dat = np.genfromtxt(args.xy,dtype='float64',skip_header=args.skip,usecols=(4,5),names="x,y",invalid_raise=False)

# Get number of lines...
sfile = open(args.sint,'r')
ns = len(sfile.readlines())
sfile.close()

# Create shuffler.
si = np.arange(0,ns)
rdm.shuffle(si)

# Do the thing.
sfile  = open(args.sint,'r')
ofile = [open(args.sint + str(i+1),'w') for i in np.arange(0,args.nsplit)]
  
#ofile1 = open(args.sint + '1','w')
#ofile2 = open(args.sint + '2','w')

for i in np.arange(0,ns):
  iw = rdm.choice(dat)
  ofile[si[i] % args.nsplit].write( "%s%8.1f%9.1f\n" % (sfile.readline().splitlines()[0], iw[0], iw[1] ))

sfile.close()
for i in np.arange(0,args.nsplit):
  ofile[i].close()
  
