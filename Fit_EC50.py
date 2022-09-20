#!/usr/bin/python
import sys
import os
import os.path
import commands
import time
import math
import copy
import pprint
import collections
import numpy as np
from collections import OrderedDict
from scipy.optimize import curve_fit
import scipy
from scipy import stats
from scipy.stats import norm
from numpy import random
#from numpy.random import Generator


try:
        infilename1   = sys.argv[1]
     #   outfilename   = sys.argv[3]
        inital_guess  = [float(sys.argv[2]), float(sys.argv[3]) ]

    #    ns_per_step  = float (sys.argv[6])
    #    try:
    #        ratio        =  float (sys.argv[7])
    #    except:
    #        ratio        =  1 
except:
     print "Usage:",sys.argv[0], "infile  ns per step  outfile"; sys.exit(1)

#########################################


print('       ')
print('       ')
print('       ')
print('       ')
print('       ')


ifile1 = open(infilename1,'r') # open  index pdb, we need resid 


ofile = open(str(infilename1+'.fit'),'w') # open file for writing


def func_hill(L, EC50 , Gmax):
    return Gmax*L / (L + EC50 ) 

counter   =  0

data = OrderedDict()

L_index  = []
GV_shift_index_current = []
SEM_index = []


for line in ifile1:

    columns = line.split()

    if len(columns) == 3  : # and line[0:4] == 'ATOM':

       L_current  = float(columns[0])
       GV_shift_current = float(columns[1])
       SEM_current = float(columns[2]) 

       L_index.append(L_current)       
       GV_shift_index_current.append(GV_shift_current)
       SEM_index.append(SEM_current)
        
  
       data[L_current] = OrderedDict()

       data[L_current]['GV_shift'] = GV_shift_current
       data[L_current]['SEM']      = SEM_current
 
print 'if not add SEM'
popt, pcov  = curve_fit(func_hill, L_index, GV_shift_index_current)

#popt array
#    Optimal values for the parameters so that the sum of the squared residuals of f
#pcov a2-D array
#    The estimated covariance of popt. The diagonals provide the variance of the parameter estimate. To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov)). 
#

print 'EC50, Gmax, fit with raw mean ', popt

print 'if  add SEM'

popt, pcov  = curve_fit(func_hill, L_index, GV_shift_index_current, inital_guess, SEM_index, True)

print 'EC50, Gmax, fit with raw mean ', popt

#print 'estimated covariance of popt', pcov

perr = np.sqrt(np.diag(pcov))

print perr

##########   

print ' now  generate a data set with same mean on each L'


data_reconstruction = OrderedDict()


for L_current in L_index:
   GV_shift_current = data[L_current]['GV_shift']
   SEM_current      = data[L_current]['SEM'] 
  
   num_samples = 1000

   data_reconstruction[L_current] = OrderedDict()
 
   data_reconstruction[L_current]['GV_shift'] =  np.random.normal(GV_shift_current, SEM_current, num_samples)
 #  data_reconstruction[L_current]['GV_shift'] =  np.random.Generator.normal(GV_shift_current, SEM_current, num_samples)    # we do not have thie model in current version


   print L_current, np.mean(data_reconstruction[L_current]['GV_shift']),  np.std(data_reconstruction[L_current]['GV_shift'])


###########   now 

print ' not fit from random choice  '

# each time  GV_shift_index_current is regenerate 

num_samples = 10000

current_samples = 0

Gmax_list = []
EC50_list = []

while current_samples < num_samples:

    GV_shift_index_current = []

    for L_current in L_index:
    
      GV_shift_current = random.choice( data_reconstruction[L_current]['GV_shift'] )
      GV_shift_index_current.append(GV_shift_current)

    popt, pcov = curve_fit(func_hill, L_index, GV_shift_index_current)
    
    Gmax_list.append(popt[1])
    EC50_list.append(popt[0])

    if float(current_samples)%1000 == 0 :
       print  ' done frame ', current_samples

    current_samples += 1 

print 'EC50 sampled',  np.mean(EC50_list),  np.std(EC50_list)   
print 'Gmax sampled',  np.mean(Gmax_list),  np.std(Gmax_list)













