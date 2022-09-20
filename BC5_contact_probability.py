#!//usr/bin/python

import MDAnalysis
import numpy
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
from collections import OrderedDict


structurefile = sys.argv[1]
targettrajectory = sys.argv[2]

outfilename   = sys.argv[3]

try:
    intitial_frame = int(sys.argv[4])
    final_frame    = int(sys.argv[5])   

except:
    intitial_frame = 0
    final_frame    = 99999999

contact_dist = 6

print  " "
print  " "
print  " the contact probability of BC5 to BK is calculated, then averaged by 4 chain and time "
print  " "
print  " note the chain order is important, check chain order section before use  "
print  " "
print  " "
print  " "



print  "  NOTE!  insertion / deletion need make a local version as resid changed !!  "

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)

#gateresidue = targetprotein.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")


############################   for  inertion/deletion mutations, modify above range is enough   all  folowin part call this range


##############  chaeck 207 210 213

ofile  = open(str(outfilename + '_av_and_sep.xvg'),'w') # open file for writing


RCKchain = ['B','C','D','A']   ## may be wrong ....
VSDchain = ['A','B','C','D']  

all_residue_prob = OrderedDict()
all_residue_list = []

# initalize value
framenumber = 0

framenumber_effective = 0.0

for ts in targetprotein.trajectory:
   framenumber += 1

   if ts.frame > final_frame :

      break

   if intitial_frame != 0  and ( framenumber <= intitial_frame or framenumber >= final_frame ):
      pass

   elif intitial_frame == 0 or ( framenumber > intitial_frame and framenumber < final_frame ):

      #x0 = gateresidue.center_of_mass()[0]
      #y0 = gateresidue.center_of_mass()[1]
      #z0  = gateresidue.center_of_mass()[2]


      framenumber_effective += 1

      BC5   = targetprotein.select_atoms(" resname BC5 " )
     # print BC5     
   

      for current_BC5 in  BC5.residues:
          temp_resi   = current_BC5.resid

          Current_sel   = targetprotein.select_atoms(" ( protein and not (name H*) ) and ( around  " + str(contact_dist) + " (resid %d and (resname BC5 and not (name H*) ) ) ) "%(temp_resi) )             
         
          for r in Current_sel.residues:

              current_resi = r.resid

              if current_resi not in all_residue_prob :

                 all_residue_prob[current_resi] = 0.0
                 all_residue_list.append(current_resi)

              all_residue_prob[current_resi] += 1
       
       #  Clink_331to334   = targetprotein.select_atoms(" (resid %d:%d ) and segid %s "%(Clink_331to334_id[0], Clink_331to334_id[-1], currect_VSDchain) )
       #   diff          = Current_sel.center_of_mass()[2] -z0
          # r                    = Clink_321.center_of_mass()-gateresidue.center_of_mass()
          # dist_321_filter     += numpy.linalg.norm(r)
          #  now there is question   how define the 

      if float(ts.frame)%100 == 0 :
          print  ' done frame ', framenumber


all_residue_list = sorted(all_residue_list)

for current_resi in all_residue_list:

 # print current_resi,  (all_residue_prob[current_resi] / 4)/framenumber_effective
   ofile.write( str(current_resi) + " " + str( (all_residue_prob[current_resi] / 4)/framenumber_effective ) + "\n")

                            ##  average by chain, then by frame ...






