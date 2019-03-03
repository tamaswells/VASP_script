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
        #self.lattice*=self.scaling_factor
        self.element_list=[j for j in poscar.readline().split()]
        try:
            self.element_amount=[int(j) for j in poscar.readline().split()]
        except ValueError:
            raise ValueError('VASP 5.x POSCAR is needed!')
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
            self.atomic_position*=self.scaling_factor
        poscar.close()
        return self.atomic_position

    def xyz_write(self):
        print('Now writing new vasp structures.')
        if os.path.exists("POSCAR.bak"):
            shutil.copyfile("POSCAR.bak","POSCAR.bak1")
            os.remove("POSCAR.bak")
        writen_lines=[]
        writen_lines.append(self.title)
        writen_lines.append(str(self.scaling_factor))
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
        poscar=open("POSCAR.bak",'w')
        poscar.writelines(writen_lines)
        poscar.close()        

    def energy_extract(self):
        print('Now reading vasp energies.')
        if os.path.exists("OUTCAR"):
            outcar=open("OUTCAR",'r')
            for index,line in enumerate(outcar):
                if "energy(sigma->0)" in line:
                    E0=float(line.split()[-1])
            outcar.close()
            self.energy=E0
            return E0         
        elif os.path.exists("OSZICAR"):
            oszicar=open("OSZICAR",'r')
            for index,line in enumerate(oszicar):
                if "E0=" in line:
                    E0=float(line.split()[4])
            oszicar.close()
            self.energy=E0
            return E0                         
        else:
            raise IOError('OSZICAR OR OUTCAR does not exist！')

if __name__ == "__main__":
    poscar=VASP()
    cartesian_position=poscar.xyz_read()
    print(cartesian_position)
    poscar.xyz_write()
