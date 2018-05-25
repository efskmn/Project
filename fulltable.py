#!/usr/bin/python3
#
############################# PROGRAM DESCRIPTION #####################################
#
# Create a big astropy table and a dictionary for all the stars from all the detectors.
#
############################# BLOCK IMPORTS AND STUFF #################################

import numpy as np                                     
from glob import glob                                  
from astropy.table import Table                        
import gc                                                       

#######################################################################################


myd = glob("C*.rawh5")                                                                         
myd.sort()                                                                                                        
print (myd)                                                                                      
ddict = {}                                                                                       
names = []; datas = []                                                                                              
for i in myd:                                                                                                       
    name = i                                                                              
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
                                                                             
def increment(basla,bitis,i, incr, to):                                                                             
    m = (to["Jerr"] < 0.2) & (to["Kerr"] < 0.2) & (to["Sharp"] < 1.0) & (to["Sharp"] > -1.0) & (to["CS"] < i+incr) & (to["CS"] > i)                                                                                                  
    print (i, "m", i+incr)                                                                                          
    x = to["J"]; y = to["K"]; c = to["CS"]                       
    return x[m], y[m], c[m]                                                                                         


l = []; b=[]; j=[]; je=[]; k=[]; ke=[]; sh=[]; ch=[]; cs=[]                                      
for i in (sorted(ddict.keys())):                                                                 
    print(i)                                                                                                       
    l.extend(ddict[i]['data']["L"])                                                                                
    b.extend(ddict[i]['data']["B"])                                                       
    j.extend(ddict[i]['data']["J"])                                                       
    je.extend(ddict[i]['data']["Jerr"])                                                                            
    k.extend(ddict[i]['data']["K"])                                                                                
    ke.extend(ddict[i]['data']["Kerr"])                                                                   
    sh.extend(ddict[i]['data']["Sharp"])                                                              
    ch.extend(ddict[i]['data']["Chi"])                                                          
    #cs.extend(ddict[i]['data']["CS"])                                                                              

to = Table([l,b,j,je,k,ke,sh,ch], names=["l","b","j","je","k","ke","sh","ch"])#,"cs"])