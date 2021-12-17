#!/usr/bin/env python

import numpy as np
import sys
import math
import MDAnalysis 
from MDAnalysis.analysis.distances import dist

u = MDAnalysis.Universe('XDATCAR.pdb', dt=1.0) 

donors = u.select_atoms('name O')
hydrogens = u.select_atoms('name H')
acceptors = u.select_atoms('name O')

cutoff = 3.5
angle = 30
H_per_donor=2

#http://jerkwin.github.io/2016/12/31/GROMACS%E5%92%8CVMD%E4%B8%AD%E7%9A%84%E6%B0%A2%E9%94%AE%E5%88%A4%E5%AE%9A%E6%A0%87%E5%87%86/
def smallest_distance_to(A, group_B, box):
    return dist(MDAnalysis.core.groups.AtomGroup([A for i in range(len(group_B))]), group_B, box=box)[-1]

def _angle_between(r_ij, r_ik, r_jk):
    costheta = (0.5/(r_ij*r_ik) * (r_ij * r_ij + r_ik* r_ik - r_jk * r_jk))[0]
    if costheta<-1:
        costheta=-1
    elif costheta>1:
        costheta=1
    return math.degrees(math.acos(costheta))
    
total_bonds = 0
for ts in u.trajectory[:]:
    count=0
    #print ts.frame
    matches = dict()
    for acc in acceptors:
        # distance criterion
        acc_donors_dist = smallest_distance_to(acc, donors, ts.dimensions)

        acceptable_donors, = np.where((acc_donors_dist <= cutoff) & (acc_donors_dist > 0))
        
        # angle criterion
        for candidate_donor in acceptable_donors:
            hydrogen_donors_dist = smallest_distance_to(donors[candidate_donor], hydrogens, ts.dimensions)

            acceptable_hydrogens = np.argpartition(hydrogen_donors_dist, H_per_donor-1)[:H_per_donor] # donor_hydrogen distance  
            r_D_A = acc_donors_dist[candidate_donor]
            for candidate_hydrogen in acceptable_hydrogens:
                r_D_H = hydrogen_donors_dist[candidate_hydrogen]
                r_H_A = smallest_distance_to(hydrogens[candidate_hydrogen], MDAnalysis.core.groups.AtomGroup([acc]),ts.dimensions)
                alpha = _angle_between(r_D_A, r_D_H, r_H_A)
                if alpha <= angle:    
                    count+=1
    total_bonds += count
    print("Now frame %9.4f, %d hydrogen bonds;" %(ts.frame, count))
print('Average %.1f.' %(total_bonds/len(u.trajectory)))


    