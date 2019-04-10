#!/bin/python
# -*- coding:utf-8 -*-
#Convert XDATCAR to PDB and extract energy & temperature profile for AIMD simulations
# -h for help;-b for frame started;-e for frame ended;;-p for PDB conversion  
#by nxu tamas@zju.edu.cn
#version 1.2
#date 2019.4.9

import os
import shutil
import numpy as np
import matplotlib as mpl
import math
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
from optparse import OptionParser  
import sys
from copy import deepcopy

class debug(object):
    def __init__(self,info='debug'):
        self.info=info
    def __call__(self,func):
        def wrapper(*args,**kwargs):
            print('[{info}] Now entering function {function}.....' \
                .format(info=self.info,function=getattr(func,"__name__")))
            return func(*args,**kwargs)
        return wrapper

class Energy_Temp(object):
    def __init__(self):
        'Read vasp MD energies and temperature.'
        self.temp=[]
        self.energy=[]
        self.energy_extract()
            
    def energy_extract(self):
        print('Now reading vasp MD energies and temperature.')        
        if os.path.exists("OSZICAR"):
            oszicar=open("OSZICAR",'r')
            for index,line in enumerate(oszicar):
                if "E0=" in line:
                    self.temp.append(float(line.split()[2]))
                    self.energy.append(float(line.split()[4]))
            oszicar.close()
            #return (self.temp,self.energy)                         
        else:
            raise IOError('OSZICAR does not exist!')

class XDATCAR(Energy_Temp):
    def __init__(self):
        'Read lattice information and atomic coordinate.'
        super(XDATCAR,self).__init__()
        print('Now reading vasp XDATCAR.')       
        self.lattice=np.zeros((3,3))
        self.NPT=False
        self.frames=0
        self._timestep=1  # timestep 1fs
        self.alpha=0.0
        self.beta=0.0
        self.gamma=0.0
        self.current_frame=0
        self.total_elements=[]
        self._lattice1=0.0
        self._lattice2=0.0
        self._lattice3=0.0
        self._lattice=np.zeros(3)
        self.format_trans=False

        @property 
        def timestep(self):
            return self._timestep
        @timestep.setter
        def timestep(self,new_time_step):
            self._timestep=new_time_step
            
        if os.path.exists("XDATCAR"):
            self.XDATCAR=open("XDATCAR",'r')
        else:
            raise IOError('XDATCAR does not exist!')
        title=self.XDATCAR.readline().strip()
        for index,line in enumerate(self.XDATCAR):
            if "Direct" in line:
                self.frames+=1
            elif title in line:
                self.NPT=True
        #self.frames=len(['Direct' for line in self.XDATCAR if "Direct" in line])
        print('Total frames {0}, NpT is {1}'.format(self.frames,self.NPT))
        assert len(self.energy)==self.frames, \
            'Number of XDATCAR frames does not equal to that of Energy terms in OSZICAR'   
        self.XDATCAR.seek(0)
        self.lowrange=0;self.uprange=self.frames-1
        if self.NPT == False: self.lattice_read()
        
    def step_select(self,selected_step): # 't>100ps and t < 1000ps'
        'eg. t > 100 and t < 1000'
        assert isinstance(selected_step,str),'Selected timestep must be in a "string"'
        if 'and' in selected_step:
            conditions=selected_step.split("and")
        else:
            conditions=[selected_step]
        for condition in conditions:
            try:
                if '>=' in condition:
                    ranges=int(condition.split('>=')[1].strip())
                    self.lowrange=ranges-1
                elif '<=' in condition:
                    ranges=int(condition.split('<=')[1].strip())
                    self.uprange=ranges-1
                elif '>' in condition:
                    ranges=int(condition.split('>')[1].strip())
                    self.lowrange=ranges
                elif '<' in condition:
                    ranges=int(condition.split('<')[1].strip())
                    self.uprange=ranges-2
                else:
                    print('Selected timestep is invaid!');continue
            except ValueError:
                print('Selected timestep is invaid!') 
                sys.exit(0)
            if (self.lowrange >= self.frames-1) or (self.uprange < 0):
                raise ValueError('Selected timestep is wrong!')
            if self.lowrange < 0: self.lowrange= 0
            if self.uprange > self.frames-1: self.uprange= self.frames-1          
                
    def lattice_read(self):    
        self.title=self.XDATCAR.readline().rstrip('\r\n').rstrip('\n')
        self.scaling_factor=float(self.XDATCAR.readline())
        for i in range(3):
            self.lattice[i]=np.array([float(j) for j in self.XDATCAR.readline().split()])
        self.lattice*=self.scaling_factor
        self._lattice1=np.sqrt(np.sum(np.square(self.lattice[0])))
        self._lattice2=np.sqrt(np.sum(np.square(self.lattice[1])))
        self._lattice3=np.sqrt(np.sum(np.square(self.lattice[2])))
        self._lattice[0]=self._lattice1
        self._lattice[1]=self._lattice2
        self._lattice[2]=self._lattice3
        self.alpha=math.acos(np.dot(self.lattice[1],self.lattice[2]) \
            /float((self._lattice2*self._lattice3)))/np.pi*180.0
        self.beta=math.acos(np.dot(self.lattice[0],self.lattice[2]) \
            /float((self._lattice1*self._lattice3)))/np.pi*180.0
        self.gamma=math.acos(np.dot(self.lattice[0],self.lattice[1]) \
            /float((self._lattice1*self._lattice2)))/np.pi*180.0
        self.element_list=[j for j in self.XDATCAR.readline().split()]
        try:
            self.element_amount=[int(j) for j in self.XDATCAR.readline().split()]
        except ValueError:
            raise ValueError('VASP 5.x XDATCAR is needed!')
        for i in range(len(self.element_amount)):
            self.total_elements.extend([self.element_list[i]]*self.element_amount[i])
        #print(self.total_elements)
        self.total_atom=sum(self.element_amount)
        self.atomic_position=np.zeros((self.total_atom,3))

    def __iter__(self):
        return self

    def __next__(self):
        self.next()

    def skiplines_(self):
        if self.NPT == True: self.lattice_read() 
        self.XDATCAR.readline()
        for i in range(self.total_atom):
            self.XDATCAR.readline()

    def next(self):
        if self.NPT == True: self.lattice_read() 
        self.XDATCAR.readline()
        for i in range(self.total_atom):
            line_tmp=self.XDATCAR.readline()
            self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])         
        self.cartesian_position=np.dot(self.atomic_position,self.lattice)
        self.current_frame+=1
        #
        #self.cartesian_position*=self.scaling_factor
        return self.cartesian_position

    def writepdb(self):
        tobewriten=[]
        tobewriten.append("MODEL         %r" %(self.current_frame))
        tobewriten.append("REMARK   Converted from XDATCAR file")
        tobewriten.append("REMARK   Converted using VASPKIT")
        tobewriten.append('CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{4:7.2f}{5:7.2f}' .format(self._lattice1,\
            self._lattice2,self._lattice3,self.alpha,self.beta,self.gamma))
        for i in range(len(self.total_elements)):
            tobewriten.append('%4s%7d%4s%5s%6d%4s%8.3f%8.3f%8.3f%6.2f%6.2f%12s' \
                %("ATOM",i+1,self.total_elements[i],"MOL",1,'    ',self.cartesian_position[i][0],\
                self.cartesian_position[i][1],self.cartesian_position[i][2],1.0,0.0,self.total_elements[i]))
        tobewriten.append('TER')
        tobewriten.append('ENDMDL')
        with open('XDATCAR.pdb','a+') as pdbwriter:
            tobewriten=[i+'\n' for i in tobewriten]
            pdbwriter.writelines(tobewriten)

    def unswrapPBC(self,prev_atomic_cartesian):
        diff= self.cartesian_position-prev_atomic_cartesian
        prev_atomic_cartesian=deepcopy(self.cartesian_position) 
        xx=np.where(diff[:,0]>(self._lattice1/2),diff[:,0]-self._lattice1\
            ,np.where(diff[:,0]<-(self._lattice1/2),diff[:,0]+self._lattice1\
                ,diff[:,0]))
        yy=np.where(diff[:,1]>(self._lattice2/2),diff[:,1]-self._lattice2\
            ,np.where(diff[:,1]<-(self._lattice2/2),diff[:,1]+self._lattice2\
                ,diff[:,1]))
        zz=np.where(diff[:,2]>(self._lattice3/2),diff[:,2]-self._lattice3\
            ,np.where(diff[:,2]<-(self._lattice3/2),diff[:,2]+self._lattice3\
                ,diff[:,2]))
        xx=xx.reshape(-1,1);yy=yy.reshape(-1,1);zz=zz.reshape(-1,1)
        return (prev_atomic_cartesian,np.concatenate([xx,yy,zz],axis=1))

    def reset_cartesian(self,real_atomic_cartesian,center_atom):
        if center_atom > len(real_atomic_cartesian)-1:
            raise SystemError("Selected atom does not exist!")
        for i in range(0,len(real_atomic_cartesian)):
            for j in range(3):
                if (real_atomic_cartesian[i][j]-real_atomic_cartesian[center_atom][j])>self._lattice[j]/2:
                    real_atomic_cartesian[i][j]-=self._lattice[j]
                if (real_atomic_cartesian[i][j]-real_atomic_cartesian[center_atom][j])<-self._lattice[j]/2:
                    real_atomic_cartesian[i][j]+=self._lattice[j]                
        return real_atomic_cartesian
                    
    def __call__(self,selected_step):
        self.step_select(selected_step)
    def __str__(self):
        return ('Read lattice information and atomic coordinate.')
    __repr__=__str__
    
class plot(object):
    'Plot MD temperature and energy profile'
    def __init__(self,lwd,font,dpi,figsize,XDATCAR_inst=None):
        self.XDATCAR_inst=XDATCAR_inst;self.timestep=self.XDATCAR_inst.timestep
        self.time_range=(self.XDATCAR_inst.lowrange,self.XDATCAR_inst.uprange+1)
        self.lwd=lwd;self.font=font;self.dpi=dpi;self.figsize=figsize

    #@debug(info='debug')
    def plotfigure(self,title='MD temperature and energy profile'):
        from matplotlib import pyplot as plt
        self.newtemp=self.XDATCAR_inst.temp;self.newenergy=self.XDATCAR_inst.energy
        xdata=np.arange(self.time_range[0],self.time_range[1])*self.timestep
        axe = plt.subplot(121)
        self.newtemp=self.newtemp[self.time_range[0]:self.time_range[1]]
        axe.plot(xdata,self.newtemp, \
                 color='black', lw=self.lwd, linestyle='-', alpha=1)
        with open("Temperature.dat",'w') as writer:
            writer.write("Time(fs)   Temperature(K)\n")
            for line in range(len(xdata)):
                writer.write("{0:f}  {1:f}\n" .format(xdata[line],self.newtemp[line]))
        axe.set_xlabel(r'$Time$ (fs)',fontdict=self.font)
        axe.set_ylabel(r'$Temperature$ (K)',fontdict=self.font)
        axe.set_xlim((self.time_range[0]*self.timestep, self.time_range[1]*self.timestep))
        axe.set_title('MD temperature profile')
        
        axe1 = plt.subplot(122)
        self.newenergy=self.newenergy[self.time_range[0]:self.time_range[1]]
        axe1.plot(xdata,self.newenergy, \
                 color='black', lw=self.lwd, linestyle='-', alpha=1)
        with open("Energy.dat",'w') as writer:
            writer.write("Time(fs)   Energy(ev)\n")
            for line in range(len(xdata)):
                writer.write("{0:f}  {1:f}\n" .format(xdata[line],self.newenergy[line]))        
        axe1.set_xlabel(r'$Time$ (fs)',fontdict=self.font)
        axe1.set_ylabel(r'$Energy$ (ev)',fontdict=self.font)
        axe1.set_xlim((self.time_range[0]*self.timestep, self.time_range[1]*self.timestep))
        axe1.set_title('MD energy profile')
        #plt.suptitle(title)
        fig = plt.gcf()
        fig.set_size_inches(self.figsize)
        plt.tight_layout() 
        plt.savefig('ENERGY.png',bbox_inches='tight',pad_inches=0.1,dpi=self.dpi)
        
        
if __name__ == "__main__":
    parser = OptionParser()  
    parser.add_option("-b", "--begin",  
                       dest="begin", default=1,  
                      help="frames begin") 

    parser.add_option("-e", "--end",  
                      dest="end", default='false',  
                      help="frames end") 

    parser.add_option("-t", "--timestep",  
                      dest="timestep", default=1.0,  
                      help="timestep per frame") 

    parser.add_option("-p", "--pdb",  
                      action="store_true", dest="format_trans", default=False,  
                      help="choose whether to convert XDATCAR to PDB!") 

    parser.add_option("--pbc",  
                      action="store_true", dest="periodic", default=False,  
                      help="choose whether to swrap PBC images!") 

    parser.add_option("-i", "--index", 
                      dest="index", default=-1,  
                      help="choose which atom to center whole molecule!") 

    (options,args) = parser.parse_args()
    XDATCAR_inst=XDATCAR()
    XDATCAR_iter=iter(XDATCAR_inst)
    if options.format_trans:
        XDATCAR_inst.format_trans=True
        if os.path.exists("XDATCAR.pdb"):
            shutil.copyfile("XDATCAR.pdb","XDATCAR-bak.pdb")
            os.remove("XDATCAR.pdb")
    else:
        XDATCAR_inst.format_trans=False
    XDATCAR_inst.timestep=float(options.timestep)   #timestep 1fs
    if options.end == 'false':
        XDATCAR_inst('t>= %r' %(int(options.begin)))
    else:
        XDATCAR_inst('t>= %r and t <= %r' %(int(options.begin),int(options.end))) # frame 10~300  corresponding to 20~600fs
    for i in range(XDATCAR_inst.uprange+1): 
        if (i>=XDATCAR_inst.lowrange):
            cartesian_position=XDATCAR_iter.next()
            if options.format_trans == True: 
                if options.periodic == True: 
                    if i == XDATCAR_inst.lowrange:
                        real_atomic_cartesian=deepcopy(cartesian_position)
                        if int(options.index) != -1:
                            real_atomic_cartesian=XDATCAR_inst.reset_cartesian(real_atomic_cartesian,int(options.index)-1)
                        XDATCAR_inst.cartesian_position=real_atomic_cartesian
                        prev_atomic_cartesian=deepcopy(cartesian_position)
                    else:
                        prev_atomic_cartesian,diffs=XDATCAR_inst.unswrapPBC(prev_atomic_cartesian)
                        real_atomic_cartesian+=diffs
                        XDATCAR_inst.cartesian_position=real_atomic_cartesian              
                XDATCAR_inst.writepdb()
        else:
            XDATCAR_iter.skiplines_()

    newtemp=XDATCAR_inst.temp
    newenergy=XDATCAR_inst.energy
    print('Finish reading XDATCAR.')

    timestep=XDATCAR_inst.timestep
    print('Selected time-range:{0}~{1}fs'.format((XDATCAR_inst.lowrange)*timestep,\
                        (XDATCAR_inst.uprange)*timestep))
    XDATCAR_inst.XDATCAR.close()

    lwd = 0.2  # Control width of line
    dpi=300          # figure_resolution
    figsize=(5,4)   #figure_inches
    font = {'family' : 'arial', 
        'color'  : 'black',
        'weight' : 'normal',
        'size' : 13.0,
        }       #FONT_setup
    plots=plot(lwd,font,dpi,figsize,XDATCAR_inst)
    plots.plotfigure()

