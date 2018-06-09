#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Convert direc coordiation to cartesian Writen By Qiang 
#modified by nxu
#version 1.1
#时间 2018.4.25
#新增功能：自动根据所有原子的Z坐标，将原子分成若干层（差值小于1.5A），再提醒用户选择需要固定的原子层数（从底部原子开始固定）。
#对于本身为笛卡尔坐标的POSCAR ，转换功能无法使用，但是限制原子功能仍可以使用。
#用法：python dire2cart.py (自动完成POSCAR>CONTCAR的从分数坐标转为笛卡尔坐标)  或者 python dire2cart.py your_POSCAR
#新增用法: python dire2cart.py your_POSCAR any_character(任何字符) ，提醒用户选择需要固定的原子层数，完成原子的限制，默认转化为笛卡尔坐标
import sys
import os

print """
###################################
#                                 #
#for VASP 5.2 or higher versions  #
#                                 #
###################################
"""

yuzhi = 1.5 #差值小于1.5A

if len(sys.argv) <= 1:
    print(u"默认将POSCAR，CONTCAR由分数坐标转为笛卡尔坐标")
    if not os.path.isfile("POSCAR"):
        if not os.path.isfile("CONTCAR"):
            print(u"Error,POSCAR，CONTCAR皆不存在！")
            exit()
        else:
            script = sys.argv[0]
            file_to_be_converted = "CONTCAR"
            print u"正在完成CONTCAR的转换........"

    else:
        script = sys.argv[0]
        file_to_be_converted = "POSCAR"
        print u"正在完成CONTCAR的转换........"
    
else:
    if len(sys.argv) == 2:
        print(u"正在将你的文件%s由分数坐标转为笛卡尔坐标" %sys.argv[1]) 
    else:
        print(u"正在将你的文件%s由分数坐标转为笛卡尔坐标,接下来请输入你要固定的原子层数！" %sys.argv[1]) 
    script, file_to_be_converted = sys.argv[:2]



file_read = open(file_to_be_converted, 'r')

line = file_read.readlines()
a1 = float(line[2].split()[0])
a2 = float(line[3].split()[0])
a3 = float(line[4].split()[0])
b1 = float(line[2].split()[1])
b2 = float(line[3].split()[1])
b3 = float(line[4].split()[1])
z1 = float(line[2].split()[2])
z2 = float(line[3].split()[2])
z3 = float(line[4].split()[2])
if len(sys.argv) >=3 and line[8][0]  == 'c' or line[8][0]  == 'C' or line[7][0]  == 'C' or line[7][0]  == 'c':
    a2 = a3 = b1 = b3 = z1 = z2  = 0
    a1 = b2 = z3 = 1

num_atoms = sum([int(x) for x in line[6].split()])

x_cartesian = []
y_cartesian = []
z_cartesian = []
tf = []

start_num = 9 # Default: With Selected T T T, coordination starts from line 9


def determinelayers(z_cartesian):
    seq = sorted(z_cartesian)
    min = seq[0]
    layerscount = {}
    sets = [min]
    for j in range(len(seq)):
        if abs(seq[j]-min) >= yuzhi:
            min = seq[j]
            sets.append(min)

    for i in range(1,len(sets)+1):
        layerscount[i] = []            
        for k in range(len(z_cartesian)):   
            if abs(z_cartesian[k]-sets[i-1]) <= yuzhi:
                layerscount[i].append(k)

    return layerscount


def convert():
    for i in range(start_num, num_atoms + start_num):
        x_cartesian.append(float(line[i].split()[0]) * a1 + float(line[i].split()[1]) * a2 + float(line[i].split()[2]) * a3)
        y_cartesian.append(float(line[i].split()[0]) * b1 + float(line[i].split()[1]) * b2 + float(line[i].split()[2]) * b3)
        z_cartesian.append(float(line[i].split()[0]) * z1 + float(line[i].split()[1]) * z2 + float(line[i].split()[2]) * z3)
         
        if len(line[i].split()) > 3:   # if  T T T exist, there are more than 3 elements in the list line[i].split()
            tf.append((line[i].split()[3]))
        else:
            tf.append(' ')   # if there is no T T T, use space instead. 
    layerscount =determinelayers(z_cartesian)
    #print layerscount
    if len(sys.argv) >=3:
        fixedlayer = int(input("Found %d layers, choose how many layers to be fixed." %len(layerscount)))
        for i in range(1,len(layerscount)+1):
            if i <= fixedlayer: 
                for j in layerscount[i]:
                    tf[j] = " F "
            else:
                for k in layerscount[i]:
                    tf[k] = " T "
    file_out = open(file_to_be_converted+'_C', 'w')
    
    for i in range(0,7):
        file_out.write(line[i].rstrip() + '\n')  # first 7 lines are kept the same 
    if len(sys.argv) >=3 and not 'S' in line[7]:
        file_out.write('Selective\n')    
    if 'S' in line[7]:
        file_out.write(line[7].rstrip()+ '\n')  # if  T T T exists, write the Selective line 
    file_out.write('Cartesian' + '\n')          # Coordination system is Cartesian now. 

    for i in range(0,len(x_cartesian)):
        file_out.write("\t%+-3.10f   %+-3.10f   %+-3.10f   %s %s %s\n"  
        %(x_cartesian[i], y_cartesian[i], z_cartesian[i], tf[i], tf[i], tf[i]))
    
    file_out.close()
    print '-----------------------------------------------------\n'
    print 'POSCAR with Cartesian Coordiations is named as %s_C\n' %(file_to_be_converted)
    print '-----------------------------------------------------\n'

if line[7][0]  == 'S' or line[7][0]  == 's':  # # With Selected T T T, coordination starts from line 9
    start_num = 9

    if  line[8][0]  == 'D' or line[8][0]  == 'd':
        print """
This POSCAR has Direct Coordinations, Conversion is starting....

              """
        convert()
    
    elif  line[8][0]  == 'C' or line[8][0]  == 'c':
        if  len(sys.argv) >=3:    
            convert()
            print u"发现是笛卡尔坐标，只用来限制原子!"
        else:
            print "This POSCAR has Cartesian Coordinations! Process is aborted!"



else : 
    print """
----------------------------------------------------
Pay Attetion! There is no TTT in coordinations part!
----------------------------------------------------
"""
    
    start_num = 8 # without Selected, No  T T T , coordination starts from line 8 
    
    if  line[7][0]  == 'D' or line[7][0]  == 'd':
        print """
This POSCAR has Direct Coordinations, Contersion starts....

"""
        convert()
    
    elif  line[7][0]  == 'C' or line[7][0]  == 'c':
        if  len(sys.argv) >=3:
            convert()
            print u"发现是笛卡尔坐标，只用来限制原子!"
        else:
            print "This POSCAR has Cartesian Coordinations! Process is aborted!"

file_read.close()
