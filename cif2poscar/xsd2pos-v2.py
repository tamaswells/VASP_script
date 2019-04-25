#! -*-coding:utf-8 -*-
#reference:http://muchong.com/t-1663483-1-authorid-934744
#modified by nxu
#version:1.1
#updated date:2019.4.16  
#usage:python xsd2pos.py A.xsd or python xsd2pos.py to convert A.xsd or all xsds.
#tips when using this script: 
	#1.unsymmetry your model in Material Studio (Build->Symmetry->Make P1)
	#2.if your want to fix some atoms, select these atoms, Modify->Constraints
	#3.save your model document as *.xsd file.


import re
import sys
import os
from glob import glob
import shutil
import math

class VASP(object):
    def __init__(self):
        self.lattice=[]
        self.atominfo=[]
        self.Selective_infomation=False
        self.Direct_mode=True
        self.element=[]
        self.xyzs=[]
        self.fixed_fraction=False
        self.fixed_cartesian=False
        self.atominMS=[]
        self.spacegroup=False

    def dot_product(self,array1,array2):
        assert len(array1)==len(array2)
        ret=0
        for i in range(len(array2)):
            ret+=array1[i]*array2[i]
        return ret

    def xyz_read(self,filename):
        #print('Now reading from xsd file.')
        with open (filename,'r') as reader:
            for index,line in enumerate(reader):
                if "AVector" in line or "BVector" in line or "CVector" in line:
                    pattern=u'Vector=\"(.*?)\"'
                    try:
                        for i in re.findall(pattern,line,re.S):
                            self.lattice.append([float(j) for j in i.split(',')])
                    except:
                        raise SystemError('No lattice information contained!')
                elif "Components" in line:
                    self.atominfo.append(line)
                if "RestrictedProperties" in line:
                    self.fixed_fraction=True
                if "AlongAxisSlider" in line:
                    self.fixed_cartesian=True
                if "\"P1\"" in line:
                    self.spacegroup=True
            if self.fixed_cartesian==True:
                reader.seek(0)
                alllines=reader.read()
        #lattice defined in MS is different from others;  
        if self.spacegroup == False:
            raise SystemError('No support for system with symmetry!')
        lattice1=math.pow((self.lattice[0][0]**2+self.lattice[0][1]**2+self.lattice[0][2]**2),0.5)
        lattice2=math.pow((self.lattice[1][0]**2+self.lattice[1][1]**2+self.lattice[1][2]**2),0.5)
        lattice3=math.pow((self.lattice[2][0]**2+self.lattice[2][1]**2+self.lattice[2][2]**2),0.5)

        alpha_angle=math.acos(self.dot_product(self.lattice[1],self.lattice[2])/float((lattice2*lattice3)))
        beta_angle=math.acos(self.dot_product(self.lattice[0],self.lattice[2])/float((lattice1*lattice3)))
        gamma_angle=math.acos(self.dot_product(self.lattice[0],self.lattice[1])/float((lattice1*lattice2)))


        bc2 = lattice2**2 + lattice3**2 - 2*lattice2*lattice3*math.cos(alpha_angle)
        h1 = lattice1
        h2 = lattice2 * math.cos(gamma_angle)
        h3 = lattice2 * math.sin(gamma_angle)
        h4 = lattice3 * math.cos(beta_angle)
        h5 = ((h2 - h4)**2 + h3**2 + lattice3**2 - h4**2 - bc2)/(2 * h3)
        h6 = math.sqrt(lattice3**2 - h4**2 - h5**2)
        self.lattice = [[h1, 0., 0.], [h2, h3, 0.], [h4, h5, h6]]
        pattern=u'ID=\"(\d+)\"'
        for line in self.atominfo:
            if re.findall(pattern,line,re.S):
                self.atominMS.append(int(re.findall(pattern,line,re.S)[0])) 
        
        pattern=u'Components=\"(.*?)\"'
        self.elements=[re.findall(pattern,line,re.S)[0] for line in self.atominfo]

        self.element_list = list(set(self.elements))
        self.element_list.sort(key = self.elements.index)
        self.element_amount=[self.elements.count(i) for i in self.element_list]
        pattern=u'XYZ=\"(.*?)\"'
        #self._cartesian_fixed=[[0,0,0] for i in range(len(self._atomic_position))]
        self._atomic_position=[re.findall(pattern,line,re.S)[0] for line in self.atominfo]
        self._cartesian_fixed=[['F','F','F'] if "RestrictedProperties" in line else ['T','T','T']  for line in self.atominfo ]
        self._atomic_position=[[float(j) for j in i.split(',')] for i in self._atomic_position]


        if self.fixed_cartesian==True:
            pattern =u'<AlongAxisSlider.*?RestrictedProperties.*?Objects=\"(\d+)\".*?\"LocalAxis(\d+)\".*?</AlongAxisSlider>'
            find_ret=re.findall(pattern,alllines,re.S)
            #print(find_ret)
            if find_ret:
                for i in find_ret:
                    try:
                        self._cartesian_fixed[self.atominMS.index(int(i[0]))][int(i[1])-1]='F'
                    except:
                        pass # images constrain
                    
        #print(self._cartesian_fixed)
        newsorted=sorted(zip(self.elements,self._cartesian_fixed,self._atomic_position),key=lambda tmp:self.elements.index(list(tmp)[0]))
        _,self.cartesian_fixed,self.atomic_position=tuple(zip(*newsorted))

        #del self._atomic_position,self._cartesian_fixed,newsorted
        self.title="Converted by xsd2pos.py"
        self.scaling_factor=1.0
        if self.fixed_fraction==True or self.fixed_cartesian==True:
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
                suffix=" ".join(self.cartesian_fixed[i])
                tobewriten= "{0:>15.8f}{1:>15.8f}{2:>15.8f}   "+suffix
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
        #print("%r has been writen!"%(filename))       


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


