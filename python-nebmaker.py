#!/usr/bin/python
#edit by lipai@mail.ustc.edu.cn 
#reference:http://blog.sina.com.cn/s/blog_b364ab230102wj0g.html
#2016/06/22
import numpy as np
import os 
import time
step_init=0.0001
images=input("input num of images: ")
ininame=raw_input("ini structure: ")
finname=raw_input("fin structure: ")

#read ini structure
fileopen=open(ininame,'r')
ini_data=fileopen.readlines()

head=ini_data[:9]
atom_num=sum(map(int,head[6].split()))

ini_data=ini_data[9:9+atom_num]
tmp=[]
fix=[]
for i in range(atom_num):
    tmp.append(map(float,ini_data[i].split()[0:3]))
    if(ini_data[i].split()[4]=='F'): fix.append(1)
    else: fix.append(0)
pos_a=np.array(tmp)

#read fin stucture
fileopen=open(finname,'r')
fin_data=fileopen.readlines()
fin_data=fin_data[9:9+atom_num]
tmp=[]
for i in fin_data:
    tmp.append(map(float,i.split()[0:3]))
pos_b=np.array(tmp)

#correction of periodic boundary condition 
for i in range(atom_num):
        for j in range(3):
                if(pos_a[i][j]-pos_b[i][j]>0.5): pos_a[i][j]-=1
                if(pos_a[i][j]-pos_b[i][j]<-0.5): pos_b[i][j]-=1
        
#get dist matrix and linear interpolation
dist_a=np.zeros([atom_num,atom_num])
dist_b=np.zeros([atom_num,atom_num])
for i in range(atom_num):
        for j in range(atom_num):
                tmp_a=0
                tmp_b=0
                for k in range(3):
                        tmp_a+=(pos_a[i][k]-pos_a[j][k])**2
                        tmp_b+=(pos_b[i][k]-pos_b[j][k])**2
                dist_a[i,j]=np.sqrt(tmp_a)
                dist_b[i,j]=np.sqrt(tmp_b)
dist_im=np.zeros([images,atom_num,atom_num])
pos_im=np.zeros([images,atom_num,3])
for i in range(images):
        dist_im[i]=dist_a+(i+1.0)*(dist_b-dist_a)/(images+1.0)
        pos_im[i]=pos_a+(i+1.0)*(pos_b-pos_a)/(images+1.0)

#optimization using steepest descent method
pos_tmp=np.zeros([atom_num,3])
dist_tmp=np.zeros([atom_num,atom_num])
s0=np.zeros(images)
s1=np.zeros(images)
flag=np.zeros(images)
for im in range(images):
        if(flag[im]==1): continue
        step=step_init
        print "generate image "+str(im+1)
        loop=0
        while(1):
                for i in range(atom_num):  #get dist_tmp
                        for j in range(atom_num):
                                if(i==j):
                                        dist_tmp[i,j]=10
                                else:
                                        tmp=0
                                        for k in range(3):
                                                tmp+=(pos_im[im][i][k]-pos_im[im][j][k])**2
                                        dist_tmp[i,j]=np.sqrt(tmp)
                for i in range(atom_num):
                        for sigma in range(3):
                                grad=0
                                for j in range(atom_num):
                                        if(j!=i):    
                                                grad+=2*(dist_im[im][i][j]-dist_tmp[i][j])*(pos_im[im][i][sigma]-pos_im[im][j][sigma])\
                                                                *(2*dist_im[im][i][j]-dist_tmp[i][j])/dist_tmp[i,j]**5
                                pos_tmp[i][sigma]=pos_im[im][i][sigma]+step*grad
                pos_im[im]=pos_tmp
                #judge convergence
                s0[im]=s1[im]
                s1[im]=0
                for i in range(atom_num):
                        for j in range(i):
                                s1[im]+=(dist_im[im][i][j]-dist_tmp[i][j])**2/dist_tmp[i][j]**4
                loop+=1
                print "loop:%d"%loop
                if(abs(s0[im]-s1[im])<0.01):
                        print "Converged!"
                        flag[im]=1
                        break
                if(loop>1 and s1[im]>s0[im]): step=step/3

#mkdir and generate poscar files
if(images+1<10): num='0'+str(images+1)
else:       num=str(images+1)
os.system("mkdir 00")
os.system("cp %s 00/POSCAR" %ininame)
os.system("mkdir "+num)
os.system("cp "+finname+" "+num+"/POSCAR")
for i in range(images):
        if(i+1<10): num='0'+str(i+1)
        else:       num=str(i+1)
        os.system("mkdir "+num)
        data=pos_im[i].tolist()
        filename=num+"/POSCAR"
        f=file(filename,"a+")
        f.writelines(head)
        for j in range(atom_num):
                line=map(str,data[j])
                line="  ".join(line)
                if(fix[j]==1): line=line+'    F  F  F\n'
                else:          line=line+'    T  T  T\n'
                f.write(line)
        f.close()