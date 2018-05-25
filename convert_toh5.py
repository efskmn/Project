#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #####################################
#
#
# Convert file to a hdf5 table. 
#
#
############################# BLOCK IMPORTS AND STUFF #################################

import sys
import os
import numpy as np
from astropy.table import Table, Column, MaskedColumn
import pandas as pd
from glob import glob
import argparse as argp
import os, errno
import shutil

############################# BLOCK READ INPUT ########################################
parser = argp.ArgumentParser(description='Convert data to hdf5.')
parser.add_argument('-tile',help='tile: the tilename as dxxx. Example: python3 /net/nas3/popes/efsokmen/MW/SCRIPTS/convert_toh5.py -tile d037 -pathraw /media/sdc1/RDISK/d037/d037/ -pathcp /net/nas3/popes/efsokmen/MW/DATA/VVV/DISK/ ', type=str, required=True)
parser.add_argument('-pathraw',help='pathraw: path to output of DAOMASTER, preferable CALIB_LB_*.raw file. ', type=str, required=True)
parser.add_argument('-pathcp',help='pathcp: path to cp the compressed output hdf5. ', type=str, required=True)
#######################################################################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
##########################################################
args = parser.parse_args()
tilename = args.tile
pathraw = args.pathraw
pathcp = args.pathcp
os.chdir(pathraw)
print ("Working in the directory {}".format(pathraw))
##########################################################
# This to ensure that the files actually exist...
for afile in glob(pathraw + "CALIB_LB_???"+tilename+"JK.raw"):
  if not os.path.isfile(afile):
    sys.exit(bcolors.FAIL + 'File ' + afile + ' not found!' + bcolors.ENDC)
for bfile in glob(pathraw + "CALIB_LB_???"+tilename+"JK.rawh5"):
  if os.path.isfile(bfile):
    calibh5 = glob("CALIB_LB_*h5")
    calibh5.sort()
    longdir = pathcp
    try:                       
        os.makedirs(longdir+tilename)
        print ("made the directory {}".format(longdir+tilename))
    except OSError as e:              
        print ("Folder {} exists! Skipped creating the same folder.".format(longdir+tilename))
        if e.errno != errno.EEXIST:
            raise
    for i in calibh5:
      print ("Copying {} to {}".format(i, pathcp))
      shutil.copy(pathraw + i, pathcp + tilename+'/'+i)
    sys.exit("Exiting without recreating h5 file.")
    
##########################################################
def convert_toh5(data_raw):
    o = pd.read_table(data_raw,header=1,names=("ID","L","B","J","Jerr","K","Kerr","Chi","Sharp"),delim_whitespace=True,low_memory=False, skiprows=3)
    to = Table.from_pandas(o)
    #m = (~np.isnan(to["CS"]))   # to mask the NAN values
    too = to#[m]
    namefile = data_raw+"h5"
    filo = open(namefile, "w+")
    tt = Table([too["ID"],too["L"],too["B"],too["J"],too["Jerr"],too["K"],too["Kerr"],too["Chi"],too["Sharp"]], names=["ID", "L","B", "J", "Jerr", "K", "Kerr", "Chi", "Sharp"])
    tt.write(namefile, format='hdf5',path=data_raw[:-6], serialize_meta=True, compression=True,overwrite=True)
    print ("Created {}.".format(namefile))
    filo.close()    
    return 0

myraws = glob('CALIB_LB_???'+tilename+'JK.raw')
myraws.sort()

for ind,val in enumerate(myraws):
    print ("Writing {}.".format(val))
    convert_toh5(val)

#cp *.rawh5 /net/nas3/popes/efsokmen/MW/DATA/VVV/DISK/

longdir = pathcp
try:                       
    os.makedirs(longdir+tilename)
    print ("made the directory {}".format(longdir+tilename))
except OSError as e:              
    print ("Folder {} exists! Skipped creating the same folder, but copying the new files.".format(longdir+tilename))
    if e.errno != errno.EEXIST:
        raise
    
calibh5 = glob("CALIB_LB_*h5")
calibh5.sort()
for i in calibh5: 
    shutil.copy(pathraw + i, pathcp + tilename+'/'+i)
    print ("Copied {} to {}".format(pathraw + i, pathcp+tilename))

##########################################################
#def convert_toh5_cs(data_raw):
#    o = pd.read_table(data_raw,header=1,names=("ID","L","B","J","Jerr","K","Kerr","Chi","Sharp", "CS"),delim_whitespace=True,low_memory=False, skiprows=3)
#    to = Table.from_pandas(o)
#    m = (~np.isnan(to["CS"]))   # to mask the NAN values
#    too = to#[m]
#    namefile = data_raw+"h5"
#    filo = open(namefile, "w+")
#    tt = Table([too["ID"],too["L"],too["B"],too["J"],too["Jerr"],too["K"],too["Kerr"],too["Chi"],too["Sharp"],too["CS"]], names=["ID", "L","B", "J", "Jerr", "K", "Kerr", "Chi", "Sharp","CS"])
#    tt.write(namefile, format='hdf5',path=data_raw[:-6], serialize_meta=True, compression=True,overwrite=True)
#    print ("Created {}.".format(namefile))
#    filo.close()    
#    return 0
#def convert_toh5_ee(data_raw):
#    o = pd.read_table(data_raw,header=1,names=("ID","L","B","J","Jerr","K","Kerr","Chi","Sharp", "CS", "EE"),delim_whitespace=True,low_memory=False, skiprows=3)
#    to = Table.from_pandas(o)
#    m = (~np.isnan(to["CS"]))   # to mask the NAN values
#    too = to#[m]
#    namefile = data_raw+"h5"
#    filo = open(namefile, "w+")
#    tt = Table([too["ID"],too["L"],too["B"],too["J"],too["Jerr"],too["K"],too["Kerr"],too["Chi"],too["Sharp"],too["CS"], too["EE"]], names=["ID", "L","B", "J", "Jerr", "K", "Kerr", "Chi", "Sharp","CS","EE"])
#    tt.write(namefile, format='hdf5',path=data_raw[:-6], serialize_meta=True, compression=True,overwrite=True)
#    print ("Created {}.".format(namefile))
#    filo.close()    
#    return 0
#
### Read the h5:
#from astropy.table import Table
#l = []; n = glob("*h5")
#for ind, val in enumerate(n):         
#    t = Table.read(val, path=val[:17])
#    l.append(t)         
