from optparse import OptionParser
import sys
import numpy as np
import os
import shutil

def move_selected_atom(selected_atom,super_x,super_y,x_grid=5,y_grid=5):
    range_x=np.linspace(0,1.0/super_x,x_grid+1)
    range_y=np.linspace(0,1.0/super_y,y_grid+1)
    if selected_atom==-1:
        selected_atom=len(coord)
    z=coord[selected_atom-1][-1]
    for m in range_x.tolist()[:-1]:
        for n in range_y.tolist()[:-1]:
            coord[selected_atom-1][0]=m
            coord[selected_atom-1][1]=n
            selective_info[selected_atom-1][0]='F'
            selective_info[selected_atom-1][1]='F'
            selective_info[selected_atom-1][2]='T'                        
            filename="%f-%f" %(m,n)
            if not os.path.exists(filename):
                os.mkdir(filename)
            with open(filename+"/POSCAR",'w') as writer:
                writer.write("By NXU\n")
                writer.write("%f\n" %(scaling_factor))
                for i in range(3):
                    writer.write("%7.3f %7.3f %7.3f\n" %(tuple(lattice[i])))
                writer.write(' '.join(element_type)+"\n")
                writer.write(' '.join(list(map(str,element_nums)))+"\n")
                writer.write('Selective\n')
                writer.write('Direct\n')
                for i in range(coord.shape[0]):
                    writer.write("%7.3f %7.3f %7.3f %s %s %s\n" %(coord[i][0],
                        coord[i][1],coord[i][2],selective_info[i][0],\
                        selective_info[i][1],selective_info[i][2]))                

def read_lattice(files):
    title=files.readline()
    scaling_factor=float(files.readline())
    lattice=np.zeros((3,3))
    for i in range(3):
        lattice[i]=list(map(float,files.readline().split()))
    #print(lattice)
    element_type=files.readline().split()
    element_nums=list(map(int,files.readline().split()))
    #print(element_nums)
    all_element=[]
    for i in range(len(element_type)):
        for j in range(element_nums[i]):
            all_element.append(element_type[i])
    return (scaling_factor,lattice,element_type,element_nums,all_element)
    #print(all_element)

def read_coord(files):
    tmp=files.readline()
    selective=False
    if tmp.strip().upper().startswith("S"):
        selective=True
        tmp=files.readline()
    selective_info=[['T','T','T'] for i in range(len(all_element))]
    all_coordinate=np.zeros((len(all_element),3))
    if selective:
        for i in range(len(all_element)):
            lines=files.readline()
            all_coordinate[i]=list(map(float,lines.split()[:3]))
            selective_info[i]=lines.split()[3:6]
    else:
        for i in range(len(all_element)):
            lines=files.readline()
            all_coordinate[i]=list(map(float,lines.split()[:3]))        
    if tmp.strip().upper().startswith("C"):
        all_coordinate=np.dot(all_coordinate,np.linalg.inv(lattice))
    return (selective_info,all_coordinate)

def extract(super_x,super_y,selected_atom):
    alllist=os.listdir(os.getcwd()) 
    dirlist=[] 
    filelist=[]
    each_case_list=[]
    for i in alllist:
        if not os.path.isfile(i): 
            dirlist.append(i) 
        else:
            filelist.append(i) 
    for j in dirlist: 
        dir_path=os.path.join(os.getcwd(),j) 
        OUTCAR=os.path.join(dir_path,"OUTCAR")
        CONTCAR=os.path.join(dir_path,"CONTCAR")
        with open(CONTCAR,'r') as reader:
            read_lattice(reader)
            selective_info,coord=read_coord(reader)    
            if selected_atom==-1:
                selected_atom=len(coord)
            x,y,z=tuple(coord[selected_atom-1])
        all_energy=[]
        with open(OUTCAR,'r') as reader:
            for index,line in enumerate(reader):
                if '  without entropy' in line:
                    all_energy.append(line)
        energy=float(all_energy[-1].split()[6])
        each_case_list.append([x,y,z,energy])

    all_new=[]
    for i in range(super_x):
        for j in range(super_y):
            for k in each_case_list:
                all_new.append([k[0]+(float(i)/super_x),k[1]+(float(j)/super_y),k[-2],k[-1]])

    s1=np.array(all_new)
    s2=s1[s1[:,0]==0]
    s2[:,0]=1.0
    s3=s1[s1[:,1]==0]
    s3[:,1]=1.0
    point=[1.0,1.0,s1[np.where((s1[:,0]==0) & (s1[:,1]==0))][0][-2],\
    s1[np.where((s1[:,0]==0) & (s1[:,1]==0))][0][-1]]        
    all_new.extend(s2.tolist()) 
    all_new.extend(s3.tolist()) 
    all_new.append(point)           
    with open('new_direct.dat','w') as writer:
        for i in all_new:
            writer.write(' '.join(map(str,i)))
            writer.write('\n')
    all_coord=np.array(all_new)#;print(all_new)
    cartesian=np.dot(all_coord[:,0:3],lattice)
    cartesian=cartesian.tolist()
    for i in range(len(cartesian)):
        cartesian[i].append(all_new[i][-1])
    with open('new_dat-3D.dat','w') as writer:
        for i in cartesian:
            writer.write(' '.join(map(str,i)))
            writer.write('\n') 

def submit_job(pbs_file):
    alllist=os.listdir(os.getcwd()) 
    dirlist=[] 
    filelist=[]
    for i in alllist:
        if not os.path.isfile(i): 
            dirlist.append(i) 
        else:
            filelist.append(i) 
    for j in dirlist: 
        for i in filelist: 
            if ".py" not in i and "POSCAR" not in i and "CONTCAR" not in i: 
                shutil.copy(i,os.path.join(os.getcwd(),j)) 
        os.chdir(os.path.join(os.getcwd(),j)) 
        os.system(pbs_file)
        os.chdir("..")

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-g", "--generate",  
                  action="store_true", dest="generation", default=False,  
                  help="generate grids") 
    parser.add_option("-u", "--user",  
                  action="store_true", dest="user_define", default=False,  
                  help="user-defined grid density")     
    parser.add_option("-e", "--extract",  
                  action="store_true", dest="energy", default=False,  
                  help="extract energy and data") 
    (options, args) = parser.parse_args()
    
    
    super_x=4   # supercell in x
    super_y=4   # supercell in y
    x_grid=4    # grids num in each unit 
    y_grid=4
    pbs_file="qsub vasp.pbs"

    if sys.version[0]=='2':
        input=raw_input    
    if options.user_define:
        super_x=int(input('Input how many supercells in x axis-->')) 
        super_y=int(input('Input how many supercells in y axis-->')) 
        x_grid=int(input('Input interporate how many points for each unit cell in x axis,even number are suggested-->'))                        
        y_grid=int(input('Input interporate how many points for each unit cell in y axis,even number are suggested-->'))
        pbs_file=input('command to submit jobs,e.g. qsub vasp.pbs-->')
    if os.path.exists("CONTCAR"):
        file1=open("CONTCAR",'r')
    else:
        file1=open("POSCAR",'r')    
    if (options.generation and options.energy) or (not options.generation and not options.energy):
        print("Wrong parameters! -h for help" )
        sys.exit(0)
    if options.generation:
        scaling_factor,lattice,element_type,element_nums,all_element=read_lattice(file1)
        selective_info,coord=read_coord(file1)
        atom=input('Input which atom you want to move to scan adsorption energy-->')
        try:
            atom_index=int(atom)
        except:
            print("Atom index required,1 for the first atom")
            sys.exit(0)
        move_selected_atom(atom_index,super_x,super_y,x_grid,y_grid)
        with open("grid.info",'w') as writer:
            writer.write("%d %d %d\n" %(super_x,super_x,atom_index))
            for i in range(3):
                writer.write("%f %f %f\n" %(tuple(lattice[i])))
        try:
            submit_job(pbs_file)
        except:
            print("Anto submit failed!")
            sys.exit(0)
            
    if options.energy:
        if os.path.exists("grid.info"):                
            with open("grid.info",'r') as reader:
                super_x,super_x,atom_index=tuple(map(int,reader.readline().split()))  
        else:
            atom=input('Input which atom you want to move to scan adsorption energy-->')
            try:
                atom_index=int(atom)
            except:
                print("Atom index required,1 for the first atom")
                sys.exit(0)
        scaling_factor,lattice,element_type,element_nums,all_element=read_lattice(file1)
        extract(super_x,super_y,atom_index)
        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import griddata
            plt.style.use('ggplot')
            import numpy as np
            import os

            data=np.loadtxt("new_dat-3D.dat")

            grid_density=100
            grid_density_y=grid_density
            range_x=np.linspace(0,1.0,grid_density)
            range_y=np.linspace(0,1.0,grid_density)
            direct=[]
            for n in range_y:
                for m in range_x:
                    direct.append([m,n,1.0])
            lattice=np.zeros((3,3))        
            if os.path.exists("grid.info"):                
                with open("grid.info",'r') as reader:
                    reader.readline()
                    for i in range(3):
                        lattice[i]=list(map(float,reader.readline().split()))  
            elif os.path.exists("CONTCAR"):
                with open("CONTCAR") as reader:
                    reader.readline()
                    scaling=float(reader.readline()) 
                    for i in range(3):
                        lattice[i]=list(map(float,reader.readline().split()))
                    lattice=lattice*scaling           
            cartesian=np.dot(direct,lattice)
            actual_grids=cartesian[:,0:2]
            grid_z0 = griddata(data[:,0:2], data[:,3],actual_grids, method='nearest')

            plt.contourf(cartesian[:,0].reshape(grid_density_y,-1),cartesian[:,1].reshape(grid_density_y,-1),grid_z0.reshape(grid_density_y,-1),15,cmap=plt.cm.hsv,alpha=0.75)
            plt.colorbar()
            plt.show()            
        except:
            print("Plot PES error!")
    
    file1.close()
