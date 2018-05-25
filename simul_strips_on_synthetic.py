#!/usr/bin/env python3

import numpy as np
#import matplotlib
#matplotlib.use("Qt4Agg")
from glob import glob                                            
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.path as mplPath
from astropy.table import Table           
import timeit
import gc
from matplotlib.colors import LogNorm 


def chunks(l, n):                 
    for i in range(0, len(l), n): 
        yield list(l[i:i+n])      


# cd /run/media/efsokmen/Elements/RDISK/d038

myls = glob("CALIB_*erawh5")                                              
myls.sort()                                                               
adict = {}                                                 
names = []; datas = []          
for i in myls:
    name = i[9:-9]                 
    names.append(name)             
    data = Table.read(i)           
    datas.append(data)                 
for eachname in names:             
    adict['%02s' % (eachname)]  = {}    
for ind,eachdata in enumerate(datas):    
    adict[names[ind]]['data'] = eachdata
gc.collect()                          
l = []; b=[]; j=[]; je=[]; k=[]; ke=[]; sh=[]; ch=[]; cs=[]
for i in (sorted(adict.keys())):
    print(i)
    l.extend(adict[i]['data']["L"])
    b.extend(adict[i]['data']["B"])
    j.extend(adict[i]['data']["J"])
    je.extend(adict[i]['data']["Jerr"])
    k.extend(adict[i]['data']["K"])
    ke.extend(adict[i]['data']["Kerr"]) 
    sh.extend(adict[i]['data']["Sharp"]) 
    ch.extend(adict[i]['data']["Chi"]) 
    cs.extend(adict[i]['data']["CS"]) 
to = Table([l, b,j,je, k,ke,sh,ch,cs], names=["l","b","j","je","k","ke","sh","ch","cs"])   
gc.collect()



fig, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False, figsize=(9,7))                            
m = (to["je"] < 0.2) & (to["ke"] < 0.2)                                   
k=to["k"][m];j=to["j"][m]                                                                                       
cs = to["cs"][m]                                                                                                
J = j  -3.02*cs                                     
K = k  - cs                                                                             
mJ, mK = (~np.isnan(J)), (~np.isnan(K))                                                                 
h1,xedg1,yedg1,im1 = axs[0].hist2d(j-k, k, bins=(1200,1000), cmap='magma', norm=LogNorm())
h,xedg,yedg,im3 = axs[1].hist2d(J[mJ]-K[mK], K[mK], bins=(1200,1000), cmap='magma', norm=LogNorm())
#axs[0].plot(np.linspace((j[100]-k[100]), (J[100]-K[100]), 100 ), np.linspace(k[100], K[100],100), 'b-' )
#axs[1].plot(np.linspace((j[100]-k[100]), (J[100]-K[100]), 100 ), np.linspace(k[100], K[100],100), 'b-' )       
axs[0].axis([-0.2,3.8,19.13,10.62])                                                                        
axs[1].axis([-2.8,1.2,18.5,10])                                                                                 
fig.tight_layout(h_pad=0.3, w_pad=0.6, rect=(0.04,0.04,0.98,0.98))                         
fig.text(0.5,0.012, "J-Ks", ha='center',fontsize=18)                          
fig.text(0.012, 0.5, "Ks", va='center', rotation='vertical', fontsize=18) 
fig.subplots_adjust(right=0.85)                        
cbar_ax = fig.add_axes([0.935, 0.05, 0.014, 0.88])                                    
fig.colorbar(im3, cax=cbar_ax)                                                                           
ROI_3 = roipoly(roicolor='r')                          
#axs[0].plot(rr_2[0], rr_2[1], 'k-')     
#axs[1].plot(rr_3[0],rr_3[1],'k-')
                                                  

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
                                                                  
sj = tdict['s44']['data'][:,17]
sk = tdict['s44']['data'][:,19]


fig, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False, figsize=(10,7))
J = j  -3.02*cs                                                                      
K = k  - cs                                                                             
gc.collect()                                                                                                               
mJ, mK = (~np.isnan(J)), (~np.isnan(K))
mask_3 = roi(ROI_3, j-k, k)                                                             
mask_4 = roi(ROI_4, J-K, K)                                                                             
binsize = np.arange(10, 17, 0.07)
h1,xedg1,yedg1,im1 = axs[0].hist2d(j-k, k, bins=(1200,1000), cmap='magma', norm=LogNorm())
h,xedg,yedg,im3 = axs[1].hist2d(J[mJ]-K[mK], K[mK], bins=(1200,1000), cmap='magma', norm=LogNorm())
#n1,bins1,patches1 = axs[0].hist(to["k"][m][np.where(mask_3)], bins=binsize, histtype="step", lw=1.2, color='black')
#n2,bins2,patches2 = axs[1].hist(K[np.where(mask_4)], bins=binsize, histtype="step", lw=1.2, color='black')
fig.tight_layout(h_pad=0.3, w_pad=0.6, rect=(0.04,0.04,0.98,0.98))                                              
fig.text(0.5,0.012, "Ks", ha='center',fontsize=18)
fig.text(0.012, 0.5, "Counts", va='center', rotation='vertical', fontsize=18) 
axs[0].plot(ROI_3.allxpoints, ROI_3.allypoints, 'k-')
axs[1].plot(ROI_4.allxpoints, ROI_4.allypoints, 'k-')


# For the CDF we need to sort the input, in this case, color
sorted_jks = sorted(to["j"][m][np.where(mask_3)]-to["k"][m][np.where(mask_3)])
binsize1 = np.arange(0.5,3.9, 0.1)
chunked_mask_ms = list(chunks(sorted_jks-np.ones(len(sorted_jks))*1.05,len(sj)))

fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, figsize=(16,5))
axs[0].step(sorted_jks, np.arange(to["k"][m][np.where(mask_3)].size)/479763.,'ko',ms=.1)
xedges1,yedges1,patches1 =axs[1].hist(to["j"][m][np.where(mask_3)]-to["k"][m][np.where(mask_3)],bins=binsize1, histtype="step", lw=1.5, color='black',label='MS')
axs[1].step(yedges1[:-1]-1.05, xedges1, lw=1.2, color='r',label='MS shifted (1.05 mag)')
axs[0].step(np.array(sorted_jks)-1.05, np.arange(to["k"][m][np.where(mask_3)].size)/479763.,'r--',lw=1.)
axs[1].legend()
axs[1].set_xlabel("J-Ks",fontsize=13)
axs[0].set_xlabel("J-Ks",fontsize=13)
axs[0].set_ylabel("Cumulative Density Function",fontsize=13)
axs[1].set_ylabel("Counts",fontsize=13)
fig.text(0.5,0.95, "MS strip from d038 (unreddened) Right panel=1 Gyr,[Fe/H]= 0.0 dex",ha='center', fontsize=14)



for i in chunked_mask_ms[:-1]:
    random.shuffle(i)
    gc.collect()
    redderj = sj + 3.02*np.ones(len(i))*i
    redderk = sk + i
    redderj = sj + 3.02*np.ones(len(i))*i                 
    redderk = sk + i                                   
    axs[2].plot(redderj-redderk, redderk, 'k.',ms=0.03)   
axs[2].axis([-0.5,6,7,-3]) 
#axs[2].plot(redderj-redderk, redderk, 'k.',ms=0.03)
#axs[2].axis([-0.5,6,7,-3])


sj1 = tdict['s31']['data'][:,17]
sk1 = tdict['s31']['data'][:,19]
sj2 = tdict['s32']['data'][:,17]
sk2 = tdict['s32']['data'][:,19]

gc.collect()
fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, figsize=(14,7))
fig.text(0.5,0.95, "CDF of the MS strip from d038 (unreddened) applied on synthetic CMDs",ha='center', fontsize=12)
for i in chunked_mask_ms[:-1]:
    random.shuffle(i)
    redderj = sj + 3.02*np.ones(len(i))*i
    redderk = sk+ i
    redderj1 = sj1 + 3.02*np.ones(len(i))*i
    redderk1 = sk1 + i
    redderj2 = sj2 + 3.02*np.ones(len(i))*i
    redderk2 = sk2 + i
    axs[0].hist2d(redderj-redderk,redderk, bins=(600,500), cmap='magma', norm=LogNorm())
    axs[1].hist2d(redderj1-redderk1, redderk1, bins=(500,400), cmap='magma', norm=LogNorm())
    h,xe1,ye1,im3 = axs[2].hist2d(redderj2-redderk2, redderk2, bins=(500,400), cmap='magma', norm=LogNorm())
axs[0].axis([-0.6,4,6,-4])
axs[1].axis([-0.6,4,6,-4])
axs[2].axis([-0.6,4,6,-4])
axs[0].set_title("A=1 Gyr [Fe/H]=0")
axs[1].set_title("A=10 Gyr [Fe/H]=-0.2")
axs[2].set_title("A=6 Gyr [Fe/H]=-0.2")
fig.tight_layout(h_pad=0.3, w_pad=0.6, rect=(0.04,0.04,0.98,0.95))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.935, 0.05, 0.014, 0.88])
fig.colorbar(im3, cax=cbar_ax)
axs[0].set_ylabel("$M_{Ks}$",fontsize=16)
fig.text(0.5, 0.04, '$M_{J}-M_{Ks}$', ha='center', fontsize=16)
#fig.savefig("d038_MS_CDF_2Dhist_3sCMD.png",format='png')


# Simulate a gaussian distance at 8kpc with sigma of 4 kpc:

a = np.random.normal(loc = 8000, scale = 4000, size = len(sj))
xedges,yedges,pathces=axa.hist(a,bins=30)
fig, axs = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(7,7))
axs.hist(a,bins=30,histtype='step',lw=3,color='b',label="$\mu=8kpc$")
axs.axis([-6000,25500,0,2500])
axs.plot(np.zeros(len(xedges)), np.linspace(0,2500,len(xedges)), 'r--',lw=2)
axs.plot(np.ones(len(xedges))*a.std(), np.linspace(0,2500,len(xedges)),'g:',lw=2,label='$\sigma=4kpc$')
axs.plot(np.ones(len(xedges))*3*a.std(), np.linspace(0,2500,len(xedges)),'g:',lw=2)
axs.legend()
axs.set_ylabel("Counts",fontsize=15)
axs.set_xlabel("d (pc)",fontsize=15)
#fig.savefig("distance_sim_8kpc.png",format='png')

sorted_dist = sorted(a)

xedges,yedges,pathces=axa.hist(a,bins=30)
fig, axs = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False, figsize=(7,4))
axs[0].hist(a,bins=30,histtype='step',lw=3,color='b',label="$\mu=8kpc$")
axs[0].axis([-6000,25500,0,2500])
axs[0].plot(np.zeros(len(xedges)), np.linspace(0,2500,len(xedges)), 'r--',lw=2)
axs[0].plot(np.ones(len(xedges))*a.std(), np.linspace(0,2500,len(xedges)),'g:',lw=2,label='$\sigma=4kpc$')
axs[0].plot(np.ones(len(xedges))*3*a.std(), np.linspace(0,2500,len(xedges)),'g:',lw=2)
axs[0].legend()
axs[0].set_ylabel("Counts",fontsize=15)
axs[0].set_xlabel("d (pc)",fontsize=15)
axs[1].step(sorted_dist,np.arange(len(sj))/20000.,'ko',ms=.1)
axs[1].plot(np.zeros(len(xedges)), np.linspace(0,1.3,len(xedges)), 'r--',lw=2)
axs[1].axis([-6000,20000,0,1.05])
axs.set_ylabel("CDF",fontsize=15)
axs[1].set_ylabel("CDF",fontsize=15)
axs[1].set_xlabel("d (pc)",fontsize=15)
#fig.savefig("distance_sim_8kpc.png",format='png')


chunked_mask_dist = list(chunks(sorted_dist,len(sj1)))

lredderj=[];lredderk=[];lredder1j=[];lredder1k=[];lredder2j=[];lredder2k=[];
for i in chunked_mask_ms[:-1]:                            
    random.shuffle(i)               
    redderj = sj + 3.02*np.ones(len(i))*i - 5. + 5.* np.log10(chunked_mask_dist[0])
    lredderj.extend(redderj)           
    redderk = sk+ i- 5. + 5.* np.log10(chunked_mask_dist[0])
    lredderk.extend(redderk)        
    redderj1 = sj1 + 3.02*np.ones(len(i))*i - 5. + 5.* np.log10(chunked_mask_dist[0])
    lredder1j.extend(redderj1)         
    redderk1 = sk1 + i - 5. + 5.* np.log10(chunked_mask_dist[0])
    lredder1k.extend(redderk1)   
    redderj2 = sj2 + 3.02*np.ones(len(i))*i - 5. + 5.* np.log10(chunked_mask_dist[0])
    lredder2j.extend(redderj2)   
    redderk2 = sk2 + i - 5. + 5.* np.log10(chunked_mask_dist[0])
    lredder2k.extend(redderk2)       

lredderj = np.array(lredderj)
lredderk = np.array(lredderk)
lredder1j = np.array(lredder1j)
lredder1k = np.array(lredder1k)
lredder2j = np.array(lredder2j)
lredder2k = np.array(lredder2k)



fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False, figsize=(14,7))
fig.text(0.5,0.95, "CDF of the MS strip from d038 applied on synthetic CMDs, distances $\mu=8kpc, \sigma=4kpc$",ha='center', fontsize=12)   
axs[0].plot(lredderj-lredderk,lredderk,  'k.',ms=0.03)
axs[1].plot(lredder1j-lredder1k, lredder1k, 'k.',ms=0.03)
axs[2].plot(lredder2j-lredder2k, lredder2k,  'k.',ms=0.03)
axs[0].set_title("A=1 Gyr [Fe/H]=0")
axs[1].set_title("A=10 Gyr [Fe/H]=-0.2")
axs[2].set_title("A=6 Gyr [Fe/H]=-0.2")
axs[0].axis([-0.6,4,21,11])
axs[1].axis([-0.6,4,21,11])
axs[2].axis([-0.6,4,21,11])
sj_ = sj - 5. + 5.* np.log10(8000)
sk_ = sk - 5. + 5.* np.log10(8000)
sj1_=sj1- 5. + 5.* np.log10(8000)
sj2_=sj2- 5. + 5.* np.log10(8000)
sk1_=sk1- 5. + 5.* np.log10(8000)
sk2_=sk2- 5. + 5.* np.log10(8000)
axs[0].plot(sj_- sk_, sk_,'r.',ms=.5)
axs[1].plot(sj1_ - sk1_, sk1_, 'r.',ms=0.5)
axs[2].plot(sj2_-sk2_,sk2_,'r.',ms=.5)
axs[0].set_ylabel("$m_{Ks}$",fontsize=16)
fig.text(0.5, 0.04, '$m_{J}-m_{Ks}$', ha='center', fontsize=16)


lredderj = np.array(lredderj)                      
lredderk = np.array(lredderk)                                                                                                                 
lredder1k = np.array(lredder1k)                                                     
lredder2k = np.array(lredder2k)                                                                                                               
lredder1j = np.array(lredder1j)                                                                         
lredder2j = np.array(lredder2j)                                                                                                               
lredder3j = np.array(lredder3j)                                                                                                          
lredder3k = np.array(lredder3k)                                                                                                               
                                                         
mask_5 = roi(ROI_5, lredder1j-lredder1k, lredder1k)                                                                                           
mask_6 = roi(ROI_6, lredder1j-lredder1k, lredder1k)       
mask_5o = roi(ROI_5, lredder3j-lredder3k, lredder3k)                                                                                          
mask_6o = roi(ROI_6, lredder3j-lredder3k, lredder3k)
                                                                
                                                              
fig, axs = plt.subplots(nrows=1, ncols=4, sharex=False, sharey=False, figsize=(15,7))
fig.text(0.5,0.95, "CDF of the MS strip from d038 applied on synthetic CMDs, distances $\mu=8kpc, \sigma=4kpc$",ha='center', fontsize=12)
axs[0].plot(lredderj-lredderk,lredderk,  'k.',ms=0.03)
axs[1].plot(lredder1j-lredder1k, lredder1k, 'k.',ms=0.03)
axs[2].plot(lredder2j-lredder2k, lredder2k,  'k.',ms=0.03)
axs[3].plot(lredder3j-lredder3k, lredder3k,  'k.',ms=0.03)                        
                                                                                  
axs[0].set_title("A=1 Gyr [Fe/H]=0")                                              
axs[1].set_title("A=10 Gyr [Fe/H]=-0.2")                                          
axs[2].set_title("A=6 Gyr [Fe/H]=-0.2")
axs[3].set_title("A=10 Gyr [Fe/H]=-1.0")                                                                                                      
                                 
axs[0].axis([-0.6,4,21,11])      
axs[1].axis([-0.6,4,21,11])      
axs[2].axis([-0.6,4,21,11])      
axs[3].axis([-0.6,4,21,11])
                                     
sj_ = sj - 5. + 5.* np.log10(8000)         
sk_ = sk - 5. + 5.* np.log10(8000)    
sj1_=sj1- 5. + 5.* np.log10(8000)      
sj2_=sj2- 5. + 5.* np.log10(8000)
sk1_=sk1- 5. + 5.* np.log10(8000)        
sk2_=sk2- 5. + 5.* np.log10(8000)                          
sj3_=sj3- 5. + 5.* np.log10(8000)                          
sk3_=sk3- 5. + 5.* np.log10(8000)                          
                                                                  
axs[0].plot(sj_- sk_, sk_,'r.',ms=.4)                          
axs[1].plot(sj1_ - sk1_, sk1_, 'r.',ms=0.4)
axs[2].plot(sj2_-sk2_,sk2_,'r.',ms=.4)                                                                                                        
axs[3].plot(sj3_-sk3_,sk3_,'r.',ms=0.4)

axs[0].set_ylabel("$m_{Ks}$",fontsize=16)
axs[1].plot(ROI_5.allxpoints, ROI_5.allypoints, 'b-',lw=1)
axs[1].plot(ROI_6.allxpoints, ROI_6.allypoints, 'b-',lw=1)
axs[3].plot(ROI_5.allxpoints, ROI_5.allypoints, 'g-',lw=1)
axs[3].plot(ROI_6.allxpoints, ROI_6.allypoints, 'g-',lw=1) 
fig.text(0.5, 0.04, '$m_{J}-m_{Ks}$', ha='center', fontsize=16)




binsize1 = np.arange(-0.5,3.9, 0.14)
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=False, figsize=(11,8))
fig.text(0.5,0.93, "Histograms of RC and MS strips from the CDF of the MS strip from d038 applied on synthetic CMD \n Age 10 Gyr and [Fe/H]=-0.2 (blue),\n Age 10 Gyr and [Fe/H]=-1.0 (magenta) distances $\mu=8kpc, \sigma=4kpc$",ha='center', fontsize=10)
xedges,yedges,patches = axs[0][0].hist( lredder1j[np.where(mask_5)]-lredder1k[np.where(mask_5)],bins=binsize1,histtype="step",lw=1.5,label='MS strip',color="b")
xedges1,yedges1,patches1 = axs[0][1].hist( lredder1j[np.where(mask_6)]-lredder1k[np.where(mask_6)],bins=binsize1,histtype="step",lw=1.5,label='RC strip', color='b')
xedges2,yedges2,patches2 = axs[1][0].hist( lredder3j[np.where(mask_5o)]-lredder3k[np.where(mask_5o)],bins=binsize1,histtype="step",lw=1.5,label='MS strip',color="m")
xedges3,yedges3,patches3 = axs[1][1].hist( lredder3j[np.where(mask_6o)]-lredder3k[np.where(mask_6o)],bins=binsize1,histtype="step",lw=1.5,label='RC strip',color="m")
fig.tight_layout(h_pad=0., w_pad=0., rect=(0.04,0.04,0.95,0.95))
fig.text(0.5, 0.04, '$m_{J}-m_{Ks}$', ha='center',fontsize=16)
fig.text(0.012, 0.5, "Counts", ha='center',rotation='vertical', fontsize=16)
axs[0][0].legend()
axs[1][0].legend()
axs[0][1].legend()
axs[1][1].legend()
axs[0][0].grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)
axs[1][0].grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)
axs[1][1].grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)
axs[0][1].grid(b=True, which='major', color='gray', linestyle='-.', linewidth=0.5)
#fig.savefig("d038_MS_CDF_distance_4sCMD_regions_histograms.png",format='png')

