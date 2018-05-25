import matplotlib
matplotlib.use('Agg')

import sys                                            
import os
import mpi4py.MPI
from joblib import Parallel, delayed
import multiprocessing

import numpy as np
from astropy.table import Table, Column
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
from glob import glob
from astropy.table import Table, hstack, vstack
import gc
from astropy.stats import sigma_clip
from scipy.stats import gaussian_kde 

inputs = range(3)
rank = mpi4py.MPI.COMM_WORLD.Get_rank()
size = mpi4py.MPI.COMM_WORLD.Get_size()

def processInput(i):
    return i * i
num_cores = multiprocessing.cpu_count()



def norm_distr(n, sigma, error_prof):             
    r = np.random.normal(-1*sigma,sigma,n.shape)
    fit = error_prof(r)   
    isfit = r + n     
    return isfit     

def mask_above(num,numarr):        
    m = numarr[numarr[:,19] < num] 
    return m              
                      
def Dist_Redden(filt_arr, dkpc, redden_array):                
    "Gives a range of reddened filter at a given distance"
    dpc = dkpc * 1000.    
    min_nonzero = np.min(filt_arr[np.nonzero(filt_arr)])
    filt_arr[filt_arr == 0]   = min_nonzero
    lnewfilt = []                                             
    for u in redden_array:                                    
        newfilt = filt_arr - 5. + 5.* np.log10(dpc) + u       
        lnewfilt.append(newfilt)   
    return lnewfilt              

def Redden(filt_arr, dkpc, bandname):                    
    dpc = dkpc * 1000.                                   
    min_nonzero = np.min(filt_arr[np.nonzero(filt_arr)])
    filt_arr[filt_arr == 0]   = min_nonzero             
    newfilt = filt_arr - 5. + 5.* np.log10(dpc)
    return newfilt                                       

def norm(n, sigma):                                                     
    random_scale_ammounts = np.random.normal(-1*sigma,sigma,n.shape)
    offset_from_mean = random_scale_ammounts + n        
    return offset_from_mean                             
                                                          
def clip_n_fit(x, y, degree, sigmac, maxerror):                 
    xy = Table([x, y], names=["A","Aerr"])                      
    m = (xy["Aerr"] < maxerror)                         
    xy = xy[m]                                          
    clipxy = sigma_clip(xy["Aerr"], sigma=sigmac, iters=3)
    xy_m = xy[np.where(clipxy.data[~clipxy.mask])]              
    p = np.poly1d(np.polyfit(xy_m["A"], xy_m["Aerr"], degree))  
    return xy_m, p                                              
                         
def find_d(magk,magK):   # in kiloparsecs, magk : apparent mag, magK: absolute mag
   d = 10**((magk-magK+5)/5.)
   d = round(d/1000.,3)
   return d

def chunks(l, n):                 
    for i in range(0, len(l), n): 
        yield list(l[i:i+n])      
            
            
# cd /run/media/efsokmen/Elements/SYNTHETIC/METHOD_3
myls = glob("s*.dat")                          
myls.sort()           
tdict = {}                
names = []; datas = []
for i in myls:       
    name = i[:-4]                                                 
    names.append(name)                                            
    data = np.genfromtxt(i, skip_header=15, invalid_raise=False)
    datas.append(data)             
for eachname in names:           
    tdict['%02s' % (eachname)]  = {}
for ind,eachdata in enumerate(datas):
    tdict[names[ind]]['data'] = eachdata                          
gc.collect()                              

# cd /run/media/efsokmen/Elements/RDISK/d038
myd = glob("C*011*.erawh5")                 
myd.sort()                                  
print (myd)                                 
ddict = {}                                  
names = []; datas = []                      
for i in myd:                               
    name = i[9:-7]                          
    names.append(name)                      
    print (i)                      
    data = Table.read(i)         
    datas.append(data)              
for eachname in names:               
    print (eachname)                        
    ddict['%02s' % (eachname)]  = {}      
for ind,eachdata in enumerate(datas):       
    ddict[names[ind]]['data'] = eachdata    
gc.collect()                      

x = ddict["011d038JK"]["data"]["J"]; y = ddict["011d038JK"]["data"]["Jerr"]
xk = ddict["011d038JK"]["data"]["K"]; yk = ddict["011d038JK"]["data"]["Kerr"]

j_jerr = clip_n_fit(x, y, 4, 3, 0.2)
k_kerr = clip_n_fit(xk, yk, 4, 3, 0.2)

# First reddening run: reddenlist = np.arange(0.0,0.6,0.1)
# Second reddening run: reddenlist = np.arange(0.6,1.2,0.1)
def task_multipreddist(alistofdat, listofreddenings, distances):
  whichdat= alistofdat 
  #fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=(13,10))
  
  for aredden in listofreddenings:
    fig, axs = plt.subplots(nrows=2, ncols=4, sharex=True, sharey=True, figsize=(13,10))
    for ind,f in enumerate(whichdat):      
      for ine,a in enumerate(f):
        mys = str(a[1:2])
        jj = mask_above(12.,tdict[a]['data']) ; kk = mask_above(12.,tdict[a]['data'])        
        for inx,di in enumerate(distances):
          d = round(di,2)
          jnew = Dist_Redden(jj[:,17],d, aredden*2.86) ; knew = Dist_Redden(kk[:,19],d, aredden)
          base = np.max(knew)    
          for inj, jne in enumerate(jnew):           
            #print ("File={}, Distance={}, reddening={}".format(a,d,aredden[inj]))
            jnewe = norm_distr(jne, 0.2, j_jerr[1]); knewe = norm_distr(knew[inj],0.2, k_kerr[1])
            jknewe = np.vstack([jnewe,knewe])
            zjk = gaussian_kde(jknewe)(jknewe)
            axs[ind][ine].scatter(jnewe-knewe,knewe, c=zjk, s=0.25, edgecolor='',cmap=plt.get_cmap('magma'))            
            axs[ind][ine].plot(jne -knew[inj],knew[inj], 'k.', ms=0.01)
          axs[ind][ine].text(right-0.35,base,"  d={}kpc\n".format(d),horizontalalignment='left',verticalalignment='top',fontsize=8,color='black')
          axs[ind][ine].plot(np.linspace(-1.,3.,200),np.ones(200)*base, 'b-.', lw=.4)
          axs[ind][ine].text(right,bottom,"{} ".format(a),horizontalalignment='left',verticalalignment='top',fontsize=10,color='red')  
        fig.text(0.4, 0.98, "Extinction={}".format(str(aredden)[1:-1]), fontsize=14 )        
        axs[ind][ine].axis([-0.5,2.5,21.,11.])
        fig.text(0.5, 0.015, 'J-Ks', ha='center', fontsize=14)
        fig.text(0.015, 0.5, 'Ks', va='center', rotation='vertical', fontsize=14)
        
        plt.show()
        plt.tight_layout()
        plt.savefig("denemerun"+mys +"multipred"+str(aredden[0])[2:]+".png", format='png') 
        print ("Saved={}".format("denemerun"+mys +"multipred"+str(aredden[0])[2:]+".png"))
  return 0

left, width = -0.1, 2.
bottom, height = 20., 10.
right = left + width
top = bottom + height

distances = np.linspace(2,10,8,endpoint=True)

listreddening = [np.linspace(0.0,0.1,2), np.linspace(0.1,0.41,5), np.linspace(0.41, 0.9, 5)]

listof_dats = [ [['s26', 's25b','s25', 's24'],['s23b', 's23', 's22', 's21']], [['s36', 's35b','s35', 's34'],['s33b', 's33', 's32', 's31']],  [['s46', 's45b','s45', 's44'],['s43b', 's43', 's42', 's41']]]

results = Parallel(n_jobs=num_cores)(delayed(task_multipreddist)(listof_dats[i], listreddening, distances) for i in inputs)

