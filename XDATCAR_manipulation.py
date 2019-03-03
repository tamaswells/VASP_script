#!/bin/python
# -*- coding:utf-8 -*-
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt

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

            if (self.lowrange >= self.frames-1) or (self.uprange < 0):
                raise ValueError('Selected timestep is wrong!')
            if self.lowrange < 0: self.lowrange= 0
            if self.uprange > self.frames-1: self.uprange= self.frames-1          
                
    def lattice_read(self):    
        self.title=self.XDATCAR.readline().rstrip('\r\n').rstrip('\n')
        self.scaling_factor=float(self.XDATCAR.readline())
        for i in range(3):
            self.lattice[i]=np.array([float(j) for j in self.XDATCAR.readline().split()])
        #self.lattice*=self.scaling_factor
        self.element_list=[j for j in self.XDATCAR.readline().split()]
        try:
            self.element_amount=[int(j) for j in self.XDATCAR.readline().split()]
        except ValueError:
            raise ValueError('VASP 5.x XDATCAR is needed!')
        self.total_atom=sum(self.element_amount)
        self.atomic_position=np.zeros((self.total_atom,3))

    def __iter__(self):
        return self

    def __next__(self):
        self.next()

    def next(self):
        if self.NPT == True: self.lattice_read() 
        self.XDATCAR.readline()
        for i in range(self.total_atom):
            line_tmp=self.XDATCAR.readline()
            self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])         
        self.cartesian_position=np.dot(self.atomic_position,self.lattice)
        self.cartesian_position*=self.scaling_factor
        return self.cartesian_position

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

    @debug(info='debug')
    def plotfigure(self,title='MD temperature and energy profile'):
        from matplotlib import pyplot as plt
        self.newtemp=self.XDATCAR_inst.temp;self.newenergy=self.XDATCAR_inst.energy
        xdata=np.arange(self.time_range[0],self.time_range[1])*self.timestep
        axe = plt.subplot(121)
        self.newtemp=self.newtemp[self.time_range[0]:self.time_range[1]]
        axe.plot(xdata,self.newtemp, \
                 color='black', lw=self.lwd, linestyle='-', alpha=1)
        with open("Temperature.dat",'w') as writer:
            writer.write("Time(fs)   Temperature(K)")
            for line in range(len(xdata)):
                writer.write("{0:d}  {1:f}\n" .format(xdata[line],self.newtemp[line]))
        axe.set_xlabel(r'$Time$ (fs)',fontdict=self.font)
        axe.set_ylabel(r'$Temperature$ (K)',fontdict=self.font)
        axe.set_xlim((self.time_range[0], self.time_range[1]))
        axe.set_title('MD temperature profile')
        
        axe1 = plt.subplot(122)
        self.newenergy=self.newenergy[self.time_range[0]:self.time_range[1]]
        axe1.plot(xdata,self.newenergy, \
                 color='black', lw=self.lwd, linestyle='-', alpha=1)
        with open("Energy.dat",'w') as writer:
            writer.write("Time(fs)   Energy(ev)")
            for line in range(len(xdata)):
                writer.write("{0:d}  {1:f}\n" .format(xdata[line],self.newenergy[line]))        
        axe1.set_xlabel(r'$Time$ (fs)',fontdict=self.font)
        axe1.set_ylabel(r'$Energy$ (ev)',fontdict=self.font)
        axe1.set_xlim((self.time_range[0], self.time_range[1]))
        axe1.set_title('MD energy profile')
        #plt.suptitle(title)
        fig = plt.gcf()
        fig.set_size_inches(self.figsize)
        plt.tight_layout() 
        plt.savefig('ENERGY.png',bbox_inches='tight',pad_inches=0.1,dpi=self.dpi)
        
        
if __name__ == "__main__":
    XDATCAR_inst=XDATCAR()
    XDATCAR_iter=iter(XDATCAR_inst)

    XDATCAR_inst.timestep=1   #timestep 2fs
    XDATCAR_inst('t>= 0 and t < 100000') # frame 10~300  corresponding to 20~600fs
    for i in range(XDATCAR_inst.frames): 
        if (i>=XDATCAR_inst.lowrange) and (i<=XDATCAR_inst.uprange):
            cartesian_position=XDATCAR_iter.next()
        else:
            XDATCAR_iter.next()
    newtemp=XDATCAR_inst.temp
    newenergy=XDATCAR_inst.energy
    print('Finish reading XDATCAR.')

    timestep=XDATCAR_inst.timestep
    print('Selected time-range:{0}~{1}fs'.format((XDATCAR_inst.lowrange+1)*timestep,\
                        (XDATCAR_inst.uprange+1)*timestep))
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

