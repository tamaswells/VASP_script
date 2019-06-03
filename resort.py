#!/bin/python
# -*- coding:utf-8 -*-
import os
import numpy as np
import shutil

class VASP(object):
    def __init__(self):
        self.lattice=np.zeros((3,3))
        self.Selective_infomation=False
        self.Direct_mode=True

    def xyz_read(self):
        print('Now reading vasp structures.')
        if os.path.exists("POSCAR"):
            poscar=open("POSCAR",'r')
        elif os.path.exists("CONTCAR"):
            poscar=open("CONTCAR",'r')
        else:
            raise IOError('CONTCAR OR POSCAR does not exist！')

        self.title=poscar.readline().rstrip('\r\n').rstrip('\n')
        self.scaling_factor=float(poscar.readline())
        for i in range(3):
            self.lattice[i]=np.array([float(j) for j in poscar.readline().split()])
        self.lattice*=self.scaling_factor
        self.element_list=[j for j in poscar.readline().split()]
        try:
            self.element_amount=[int(j) for j in poscar.readline().split()]
        except ValueError:
            raise ValueError('VASP 5.x POSCAR is needed!')
        self.all_element=[]
        for i in range(len(self.element_list)):
            self.all_element.extend([self.element_list[i]]*self.element_amount[i] )
        #print(self.all_element)
        line_tmp=poscar.readline()
        if line_tmp.strip().upper().startswith("S"):
            self.Selective_infomation=True 
            line_tmp=poscar.readline()  
        else:# no atoms fixed
            self.Selective_infomation=False 
        if line_tmp.strip().upper().startswith("D"):
            self.Direct_mode=True
        elif line_tmp.strip().upper().startswith("C"):
            self.Direct_mode=False
        else:
            raise ValueError("POSCAR format is not correct!")
        total_atom=sum(self.element_amount)
        self.atomic_position=np.zeros((total_atom,3))
        self.Selective_TF=[]
        if self.Selective_infomation == True:   
            for i in range(total_atom):
                line_tmp=poscar.readline()
                self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])
                self.Selective_TF.append([j for j in line_tmp.split()[3:]]) 
        else:
            for i in range(total_atom):
                line_tmp=poscar.readline()
                self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])         
        if self.Direct_mode == True:
            self.atomic_position=np.dot(self.atomic_position,self.lattice)
            #self.atomic_position*=self.scaling_factor
        poscar.close()
        return self.atomic_position

    def xyz_write(self):
        print('Now writing new vasp structures.')
        if os.path.exists("POSCAR-bak"):
            shutil.copyfile("POSCAR-bak","POSCAR-bak1")
            os.remove("POSCAR-bak")
        writen_lines=[]
        writen_lines.append(self.title)
        writen_lines.append('1.0')
        for i in range(3):
            writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                .format(self.lattice[i][0],self.lattice[i][1],self.lattice[i][2]))
        writen_lines.append('  '+'  '.join(self.element_list))
        writen_lines.append('  '+'  '.join([str(j) for j in self.element_amount]))
        if self.Selective_infomation == True:
            writen_lines.append("Selective") 
        writen_lines.append("Cartesian") 
        if self.Selective_infomation == True:
            for i in range(np.size(self.atomic_position,0)):
                writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                    .format(self.atomic_position[i][0],self.atomic_position[i][1], \
                        self.atomic_position[i][2],self.Selective_TF[i][0],\
                            self.Selective_TF[i][1],self.Selective_TF[i][2]))                    
        else:
            for i in range(np.size(self.atomic_position,0)):
                writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                    .format(self.atomic_position[i][0],self.atomic_position[i][1], \
                        self.atomic_position[i][2]))
        writen_lines=[j+'\n' for j in writen_lines]
        poscar=open("POSCAR-bak",'w')
        poscar.writelines(writen_lines)
        poscar.close()    

def determinelayers(z_cartesian):
    seq = sorted(z_cartesian)
    min = seq[0]
    layerscount = {}
    sets = [min]
    for j in range(len(seq)):
        if abs(seq[j]-min) >= thresthold:
            min = seq[j]
            sets.append(min)
    for i in range(1,len(sets)+1):
        layerscount[i] = []            
        if i > 1:
            for k in range(len(z_cartesian)):   
                if abs(z_cartesian[k]-sets[i-1]) <= thresthold and not abs(z_cartesian[k]-sets[i-2]) <=thresthold:
                    layerscount[i].append(k)            
        else:
            for k in range(len(z_cartesian)):   
                if abs(z_cartesian[k]-sets[i-1]) <= thresthold:
                    layerscount[i].append(k)
    return layerscount

if __name__ == "__main__":
    thresthold=1.5 #A
    poscar=VASP()
    cartesian_position=poscar.xyz_read()
    substitute_elements=  [ 'Xe', 'La', 'Pr', 'Nd', 'Pm', 'Sm', \
   'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',\
     'Tl', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',\
      'Th', 'Pa', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', \
      'Es', 'Fm', 'Md', 'No', 'Lr']
    layerscount =determinelayers(cartesian_position[:,2].tolist())
    #print(layerscount)
    import sys
    if sys.version[0]=='2':
        input=raw_input
    fixedlayer = int(input("Found %d layers, choose how many layers want to substitute with various elements------>" %len(layerscount)))
    index=0
    new_element=[]
    new_element_index=[]
    for i in range(1,len(layerscount)+1):
        if i <= fixedlayer: 
            old_element=[poscar.all_element[m] for m in layerscount[i]]
            unique_element=list(set(old_element))
            for n in unique_element:
                new_element.append(substitute_elements[index])
                index+=1
                tmp=[]
                for index_,l in enumerate(old_element):
                    if l==n:
                        tmp.append(layerscount[i][index_])
                new_element_index.append(tmp)
        else:
            new_element.append(poscar.all_element[layerscount[i][0]])
            new_element_index.append(layerscount[i]) 
    old_atomic_position=cartesian_position.copy()
    poscar.element_list=new_element
    poscar.element_amount=[len(i) for i in new_element_index]
    index=0
    #print(new_element_index)
    for i in new_element_index:
        for j in i:
            poscar.atomic_position[index]= old_atomic_position[j]  
            index+=1            
    poscar.xyz_write()
