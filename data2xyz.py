#!/usr/bin/python
# -*- coding:utf-8 -*-
import numpy as np
import math


masses= [1.00794,4.002602,6.941,9.012182,10.811,12.011,14.00674,\
     15.9994,18.9984032,20.1797,22.989768,24.3050,26.981539,28.0855,\
     30.97362,32.066,35.4527,39.948,39.0983,40.078,44.955910,47.88,\
      50.9415,51.9961,54.93085,55.847,58.93320,58.69,63.546,65.39,\
      69.723,72.61,74.92159,78.96,79.904,83.80,85.4678,87.62,88.90585,91.224,\
      92.90638,95.94,98,101.07,102.90550,106.42,107.8682,112.411,\
      114.82,118.710,121.75,127.60,126.90447,131.29,132.90543,137.327,138.9055,140.115,\
      140.90765,144.24,145,150.36,151.965,157.25,158.92534,162.50,\
      164.93032,167.26,168.93421,173.04,174.967,178.49,180.9479,\
      183.85,186.207,190.2,192.22,195.08,196.96654,200.59,204.3833,\
      207.2,208.98037,209,210.0,222,223,226.025,227.028,232.0381,\
      231.03588,238.0289,237.048,244.,243.,247.,247.,251.,252.,\
      257.,258.,259.,260.0]


element=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',\
 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',\
  'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', \
  'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', \
  'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', \
  'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', \
  'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', \
  'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', \
  'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

#convert_data2POSCAR
xy=xz=yz=0
with open("system_after_md.data",'r') as reader:
    all_content=reader.readlines()
    for index,line in enumerate(all_content):
        if "atoms" in line:
            atomnums=int(line.split()[0])
        if "atom" in line and "types" in line:
            typesnums=int(line.split()[0])
        if "xlo" in line and "xhi" in line:
            x=float(line.split()[1])-float(line.split()[0])
        if "ylo" in line and "yhi" in line:
            y=float(line.split()[1])-float(line.split()[0])
        if "zlo" in line and "zhi" in line:
            z=float(line.split()[1])-float(line.split()[0])
        if "xy" in line and "xz" in line and "yz" in line:
            xy=float(line.split()[0])
            xz=float(line.split()[1])
            yz=float(line.split()[2]) 
        if "Atoms" in line:
            atomsbegin=index
        if "Masses" in line:
            Massesbegin=index                               
    print(atomnums,typesnums,x,y,z,xy,xz,yz)
newlattice=[[x,0,0],[xy,y,0],[xz,yz,z]]
print(newlattice)
for i in range(atomsbegin+1,len(all_content)):
    if all_content[i].strip("\r\n").strip("\n")!="":
        atomsbegin=i
        break
        
for i in range(Massesbegin+1,len(all_content)):
    if all_content[i].strip("\r\n").strip("\n")!="":
        Massesbegin=i
        break        
mass_list={}
for i in all_content[Massesbegin:Massesbegin+typesnums]:
    tmp=i.strip("\r\n").strip("\n").split()
    mass_list[tmp[0]]=tmp[1]

def judge_element(mass,masses):
    for index,i in enumerate(masses):
        if abs(mass-i)<0.1:
            return index

atoms=[]
for i in all_content[atomsbegin:atomsbegin+atomnums]:
    mass=float(mass_list[i.split()[2]])
    element_=element[judge_element(mass,masses)]
    atoms.append((element_,[float(i.split()[4]),float(i.split()[5]),float(i.split()[6])])) #将原子的元素（返回的最后一个元素），原子的坐标（第二个元素）和是否被限制的信息打包为一个元祖，加到原子信息列表里面。
#print(atoms)
with open("system.xyz",'w') as writer:
    writer.write("%d\n\n" %(len(atoms)))
    for i in atoms:
        writer.write("%s %7.3f %7.3f %7.3f\n" %(i[0],i[1][0],i[1][1],i[1][2]))
