#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #####################################
#
# Read the CDF file written by the script roipoly_select_regions.py and apply the MS strip and distance
# Run it in the directory where synthetic populations exist
# Example: in a working ipython window
# run /scratch1/efsan/SCRIPTS/do_red_distance.py -cdf "MS_strip_CDF_CALIB_LB_011d038JK.txt" -cdfdir "/run/media/efsokmen/Elements/RDISK/d038" -spop "s11.dat" -spopdir "/run/media/efsokmen/Elements/SIMULATIONS/SYNTHETIC_D_EJK/" -mu 7 -sigma 2

#run /scratch1/efsan/SCRIPTS/do_red_distance.py -cdf "MS_strip_CDF_CALIB_LB_032d038.dat" -cdfdir "/run/media/efsokmen/Elements/RDISK/d038" -dat "CALIB_LB_032d038JK.raw" -spop "disperso.b11disp" -spopdir "/run/media/efsokmen/Elements/SIMULATIONS/SYNTHETIC_D_EJK/" #-mu 6 -sigma 0.1 -red_slope 0.5

 
############################# BLOCK IMPORTS AND STUFF #################################

import sys
import os
import os.path
import argparse as argp
import numpy as np
import math
import pylab as pl
import gc
import random as random    
from string_bool import str2bool

############################# BLOCK READ INPUT ########################################

parser = argp.ArgumentParser(description='Apply MS strips, distance and reddening!')

parser.add_argument('-cdf',help='cdf: Cumulative distribution function array written in a text file.', type=str,required=True)

parser.add_argument('-cdfdir',help='cdfdir: Directory where cumulative distribution function text is.',type=str,required=True)

parser.add_argument('-dat',help='dat: the ALLSTAR/DAOPHOT output to work with.',required=True)

parser.add_argument('-spop',help='spop: name of the IAC-star output file.',required=True,type=str)
  
parser.add_argument('-spopdir',help='spopdir:Directory where the IAC-star output file is.',required=False, type=str, default = '/run/media/efsokmen/Elements/SIMULATIONS/SYNTHETIC_D_EJK/')

parser.add_argument('-mu',help='mu: center of the Gaussian for the distance dispersion (in Kpc)', required=True, type=float)
  
parser.add_argument('-sigma',help='sigma: sigma of the Gaussian for the distance dispersion (in Kpc)', required=True, type=float)

parser.add_argument('-red_slope',help='red_slope: reddening slope', required=False, type=float,default = 0.48)

parser.add_argument('-saveornot', help='saveornot: want to save the plot with its text file?', required=False, type=str2bool,nargs='?', default=False)


args = parser.parse_args()

#######################################################################################
def chunks(l, n):                 
    for i in range(0, len(l), n): 
        yield list(l[i:i+n])    

howlarge = int(1)


# Calculate the extinction in J from the reddening slope that is given;
AJ = (1+args.red_slope)/args.red_slope


os.chdir(args.cdfdir)
#print ("Changed directory to {}".format(args.cdfdir))

if not os.path.isfile("MS_strip_CDF_"+args.dat[:-6]+".dat"):
  print ('File ' + "MS_strip_CDF_"+args.dat[:-6]+".dat" + ' not found!')
  fil = [np.zeros(800),np.zeros(800)]
  fil = np.array(fil)
else:
  fil = np.genfromtxt("MS_strip_CDF_"+args.dat[:-6]+".dat")


os.chdir(args.spopdir)
#print ("Changed directory to {}".format(args.spopdir))
ss11 = np.genfromtxt(args.spop,skip_header=15, invalid_raise=True, usecols=(0,1)) # 17 19 (antonio) 14 16 (sebastian)
sj = ss11[:,0]
sk = ss11[:,1]

#agemetal = np.genfromtxt(args.spop,skip_header=15, invalid_raise=True, usecols=(10,11)) # antonio's output case
#age = agemetal[:,0]
#metal = agemetal[:,1]
#age = (np.array([age]*howlarge)).reshape((age.size*howlarge,))
#metal = (np.array([metal]*howlarge)).reshape((metal.size*howlarge,))


sj = (np.array([sj]*howlarge)).reshape((sj.size*howlarge,))
sk = (np.array([sk]*howlarge)).reshape((sk.size*howlarge,))

# Create a Gaussian for distance dispersion 
d_disper = np.random.normal(loc = args.mu*1000, scale = args.sigma*1000, size = len(sj))
sorted_dist = sorted(d_disper)
chunked_mask_dist = list(chunks(sorted_dist,len(sj)))

# Create a list of lists of synthetic population or the MS strip for dispersing, 
# divided in the same number of points as one another
lredderj=[]
lredderk=[]


lredden_ms = []
lchunked_sj, lchunked_sk = [], []
if len(fil[:,0]) >= len(sj):
    redden_ms = list(chunks(fil[:,0],len(sj)))
    lredden_ms.append(redden_ms)
    for list_reddening in lredden_ms[:-1]:                            
        random.shuffle(list_reddening)
        redderj = sj + AJ*np.ones(len(list_reddening))*list_reddening - 5. + 5.* math.log10(chunked_mask_dist[0])
        lredderj.extend(redderj)           
        redderk = sk + list_reddening - 5. + 5.* math.log10(chunked_mask_dist[0])
        lredderk.extend(redderk)    

elif len(sj) > len(fil[:,0]):    
    chunked_sj = list(chunks(sj,len(fil[:,0])))
    chunked_sk = list(chunks(sk, len(fil[:,0])))
    c_chunked_mask_dist = list(chunks(chunked_mask_dist[0], len(fil[:,0])))
    for index, a_sj in enumerate(chunked_sj[:-1]):
        random.shuffle(fil[:,0])
        
        try:
            dd = [np.float64(math.log10(adistance)) for adistance in c_chunked_mask_dist[index]]
            redderj = a_sj + AJ*np.ones(len(fil[:,0]))*fil[:,0] - 5. + 5.*np.ones(len(fil[:,0]))* dd
            lredderj.extend(redderj)
            redderk = chunked_sk[index] + fil[:,0] - 5. + 5.*np.ones(len(fil[:,0]))* dd
            lredderk.extend(redderk)
        except (ValueError):    
            pass       
    
else:
    pass



# Convert the lists to arrays
lredderj = np.array(lredderj)
lredderk = np.array(lredderk)

gc.collect()

fig, axs = pl.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(9,10), dpi=80, facecolor='w', edgecolor='k')
fig.text(0.5,0.93, "synthetic CMD \n distances $\mu={}kpc, \sigma={}kpc, Reddening slope = {}$".format(args.mu,args.sigma,args.red_slope),ha='center', fontsize=12)
axs.plot(lredderj-lredderk,lredderk,  'k.',ms=0.1)
axs.set_title("File: {}".format(args.spop), fontsize=14)
axs.set_xlabel("J-Ks",fontsize=16)
axs.set_ylabel("Ks",fontsize=16)
axs.set_ylim(12,20.)
axs.set_xlim(0.,2.5)
axs.invert_yaxis()
pl.grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.6)

pl.show()

mlen, mmu, msig, mslope = howlarge, int(round(args.mu,0)),int(round(args.sigma,1)),int(round(args.red_slope,2))


si = np.arange(0,len(lredderj))
random.shuffle(si)


# Save the final plot and the corresponding .dat file or not. .dat file to be used for reddening and distance dispersions
if args.saveornot is True:
    
    fig.savefig("Disp_"+str(args.spop)[:-4]+args.dat[9:-6]+"{}mu{}sigma{}slope{}".format(mlen,mmu, msig, mslope)+".png",format='png')
    print ("Saved","Disp_"+str(args.spop)[:-4]+args.dat[9:-6]+"{}mu{}sigma{}slope{}".format(mlen,mmu, msig, mslope)+".png")
    
    filow = open("Disp_"+str(args.spop)[:-4]+args.dat[9:-6]+"{}mu{}sigma{}slope{}".format(mlen,mmu, msig, mslope)+".dat", "w+")
    #filow.write("#    J-Ks    Ks        Age        Fe/H\n")
    for i in np.arange(0,len(lredderj)): 
        filow.write( "%7.3f%7.3f \n" % (lredderj[i], lredderk[i])) 
    
    filow.close()
    print ("Saved file", "Disp_"+str(args.spop)[:-4]+args.dat[9:-6]+"{}mu{}sigma{}slope{}".format(mlen,mmu, msig, mslope)+".dat")
else:
    pass


