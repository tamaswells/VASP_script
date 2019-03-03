#! -*-coding:utf-8 -*-
#reference:http://muchong.com/t-1663483-1-authorid-934744
#modified by nxu
#version:1.0
#updated date:2019.3.2  
#usage:python xsd2pos.py A.xsd or python xsd2pos.py to convert A.xsd or all xsds.
#tips when using this script: 
	#1.unsymmetry your model in Material Studio (Build->Symmetry->Make P1)
	#2.if your want to fix some atoms, select these atoms, Modify->Constraints->fix fractional position
	#3.save your model document as *.xsd file.


import re
import sys
import os
from glob import glob
import shutil


class VASP(object):
    def __init__(self):
        self.lattice=[]
        self.atominfo=[]
        self.Selective_infomation=False
        self.Direct_mode=True
        self.element=[]
        self.xyzs=[]
        self.restricted=[]
        self.atomic_position=[]

    def xyz_read(self,filename):
        #print('Now reading from xsd file.')
        with open (filename,'r') as reader:
            for index,line in enumerate(reader):
                if "Vector" in line:
                    pattern=u'Vector=\"(.*?)\"'
                    for i in re.findall(pattern,line,re.S):
                        self.lattice.append([float(j) for j in i.split(',')])
                elif "Components" in line:
                    self.atominfo.append(line)

        #lattice defined in MS is different from others;  
        #print(self.lattice)
        pattern=u'Components=\"(.*?)\"'
        self.elements=[re.findall(pattern,line,re.S)[0] for line in self.atominfo]
        self.element_list = list(set(self.elements))
        self.element_list.sort(key = self.elements.index)
        self.element_amount=[self.elements.count(i) for i in self.element_list]
        pattern=u'XYZ=\"(.*?)\"'
        self._atomic_position=[re.findall(pattern,line,re.S)[0] for line in self.atominfo]
        self._restricted=[1 if "RestrictedProperties" in line else 0  for line in self.atominfo ]
        self._atomic_position=[[float(j) for j in i.split(',')] for i in self._atomic_position]
        for i in self.element_list:
            for j in range(len(self.elements)):
                if (i==self.elements[j]):
                    self.restricted.append(self._restricted[j])
                    self.atomic_position.append(self._atomic_position[j])            
        del self._atomic_position,self._restricted
        self.title="Converted by xsd2pos.py"
        self.scaling_factor=1.0
        if self.restricted.count(1)>0:
             self.Selective_infomation=True  
        else:# no atoms fixed
             self.Selective_infomation=False 
        return 0

    def xyz_write(self,filename="POSCAR"):
        #print('Now writing vasp structure.')
        if os.path.exists("POSCAR"):
            shutil.copyfile("POSCAR","POSCAR.bak")
            os.remove("POSCAR")
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
        writen_lines.append("Direct") 
        if self.Selective_infomation == True:
            for i in range(len(self.atomic_position)):
                tobewriten= "{0:>15.8f}{1:>15.8f}{2:>15.8f}   F  F  F" \
                if self.restricted[i] == 1 else "{0:>15.8f}{1:>15.8f}{2:>15.8f}   T  T  T"
                writen_lines.append(tobewriten \
                    .format(self.atomic_position[i][0],self.atomic_position[i][1], \
                        self.atomic_position[i][2]))                    
        else:
            for i in range(len(self.atomic_position)):
                writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                    .format(self.atomic_position[i][0],self.atomic_position[i][1], \
                        self.atomic_position[i][2]))
        writen_lines=[j+'\n' for j in writen_lines]
        poscar=open(filename,'w')
        poscar.writelines(writen_lines)
        poscar.close() 
        print("%r has been writen!"%(filename))       


if __name__ == "__main__":
    if len(sys.argv)==1:
        print("Converted all xrds files by default~")
        if len(glob("*.xsd"))==0:
            raise IOError('No xrds files！Please type in filename to be converted in command line!') 
        else:
            for i in glob("*.xsd"):
                poscar=VASP()
                poscar.xyz_read(i)
                poscar.xyz_write(i.split(".")[0]+"_POSCAR")

    else:
        if os.path.exists(sys.argv[1]):
            xsdname=sys.argv[1]
            poscar=VASP()
            poscar.xyz_read(xsdname)
            poscar.xyz_write()
        else:
            raise IOError('Specified file does not exist！')    


