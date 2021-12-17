#!/usr/bin/env python

import math #引入数学库
import MDAnalysis 
from MDAnalysis.analysis.distances import dist #引入MDAnalysis计算最小映像距离的函数

u = MDAnalysis.Universe('XDATCAR.pdb', dt=1.0)  #载入轨迹，设置时间步为0.5 fs

oxygens = u.select_atoms('name O')    #选择所有的O原子
hydrogens = u.select_atoms('name H')  #选择所有的H原子

def smallest_distance_to(A, group_B, box):  #计算A原子和group B所有原子最小映像距离的函数
    return dist(MDAnalysis.core.groups.AtomGroup([A for i in range(len(group_B))]), group_B, box=box)[-1]

def angle_between_jik(r_ij, r_ik, r_jk):    #余弦定理计算 j-i-k的角度
    costheta = 0.5/(r_ij*r_ik) * (r_ij * r_ij + r_ik* r_ik - r_jk * r_jk)
    return math.degrees(math.acos(costheta))

def find_the_last_two_min(lists):  #找到O和所有的H的距离中最小的两个距离对应H的序号
    min = 0
    for i in range(len(lists)):
        if lists[i] < lists[min]:
            min = i
    if min == 0:
        second_min = 1
    else:
        second_min = min - 1
    for i in range(len(lists)):
        if (lists[i] < lists[second_min]) and (i != min):
            second_min = i
    return(min,second_min)
        
bond_ave_per_frame = []    #用于储存所有帧的键长
angle_ave_per_frame = []   #用于储存所有帧的角度
for ts in u.trajectory[10000:]:   #对10000帧之后的轨迹进行迭代
    bonds=[]  #用于储存每一帧中的所有的O-H键长
    angles = [] #用于储存每一帧中的所有的H-O-H键角
    for oxygen in oxygens: #对所有的O原子进行迭代
        dist_H_O = smallest_distance_to(oxygen, hydrogens, ts.dimensions) #计算一个O原子和所有H原子的最小映像距离
        two_H_of_O = find_the_last_two_min(dist_H_O) #找到O原子所连的两个H原子的序号
        r_ij = dist_H_O[two_H_of_O[0]]   #O-H1的键长
        r_ik = dist_H_O[two_H_of_O[1]]   #O-H2的键长
        bonds.append(r_ij)  #添加到键长集合中
        bonds.append(r_ik)  
        r_jk = smallest_distance_to(hydrogens[two_H_of_O[0]], MDAnalysis.core.groups.AtomGroup([hydrogens[two_H_of_O[1]]]),ts.dimensions)
                            #计算两个H之间的距离，第二个传入的参数一定要是原子集合
        alpha = angle_between_jik(r_ij, r_ik, r_jk) #根据余弦定理计算角度
        angles.append(alpha) #添加到键角集合中
    bond_ave_per_frame.append(sum(bonds)/len(bonds)) #每一帧的平均键长
    angle_ave_per_frame.append(sum(angles)/len(angles))  #每一帧的平均角度
print("Average bond length is %.2f." %(sum(bond_ave_per_frame)/len(bond_ave_per_frame)))  #对所有帧进行平均
print("Average angle is %.2f." %(sum(angle_ave_per_frame)/len(angle_ave_per_frame)))
       