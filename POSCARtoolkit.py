#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Convert direc coordiation to cartesian Writen By Qiang 
#modified by nxu
#version 1.2
#时间 2018.8.26
#新增功能：通过flag参数控制功能使用 -h 帮助  -i 指定输入文件  
# with no extra flags, for convert POSCAR from DIRECT TO CARTESIEN
#-f(without -s -c) 根据层数固定原子  -s -f  根据原子顺序（1-5 6 9）或者元素类型（Au etc.） 选择原子，然后固定
# -s -r (TTT OR FFF OR ANY) 根据原子顺序或者元素类型 选择原子，放开，固定或者任意 

from optparse import OptionParser  
import sys 
reload(sys) 
sys.setdefaultencoding('utf-8') 
import os
import re
import copy

parser = OptionParser()  
#oarser.add_option("-f", "--file", dest="filename",  
                  #help="write report to FILE", metavar="FILE")  
parser.add_option("-i", "--in",  
                  type = 'string',dest="files", default="POSCAR",  
                  help="file to be converted or fixed")  

parser.add_option("-f", "--fix",  
                  action="store_true", dest="fixed", default=False,  
                  help="choose whether to fix atoms. [without -s] will fix layers") 

parser.add_option("-r", "--relax",  
                  action="store_true", dest="relax", default=False,   
                  help="choose whether to relax atoms, -r XXX")

parser.add_option("-s", "--select",  
                  action="store_true", dest="selected", default=False,  
                  help="select atoms to fix or relax atoms. with -r or -f") 

parser.add_option("-y", "--thresthold",  
                  type = 'float',dest="thresthold", default=1.5,  
                  help="threstholds of distance to separate layers!") 
  
(options, args) = parser.parse_args()
#print options.selected
print """
              ###################################
              # For VASP 5.2 or higher versions #
              #         Author:Li,Q;Xu,N        #
              #           Verision 1.2          #
              ###################################
"""
if options.fixed and options.relax:
    print(u"\n    Error,-r -f does not coexist")
    exit()    
if options.files == 'POSCAR':
    if not os.path.isfile("POSCAR"):
        if not os.path.isfile("CONTCAR"):
            print(u"\n    Error,POSCAR or CONTCAR does not exist!")
            exit()
        else:
            file_to_be_converted = "CONTCAR"
            print u"\n    File to be handled is *****CONTCAR*****"
    else:
        file_to_be_converted = "POSCAR"
        print u"\n    File to be handled is *****POSCAR*****"
else:
    file_to_be_converted = options.files
    print u"\n    File to be handled is *****%s*****" %options.files



file_read = open(file_to_be_converted, 'r')

line = file_read.readlines()
file_read.close()

a1 = float(line[2].split()[0])
a2 = float(line[3].split()[0])
a3 = float(line[4].split()[0])
b1 = float(line[2].split()[1])
b2 = float(line[3].split()[1])
b3 = float(line[4].split()[1])
z1 = float(line[2].split()[2])
z2 = float(line[3].split()[2])
z3 = float(line[4].split()[2])

if  line[8].strip().upper().startswith('C') or line[7].strip().upper().startswith('C') and (options.fixed or options.relax):
    a2 = a3 = b1 = b3 = z1 = z2  = 0
    a1 = b2 = z3 = 1


num_atoms = sum([int(x) for x in line[6].split()])

x_cartesian = []
y_cartesian = []
z_cartesian = []


start_num = 9 # Default: With Selected T T T, coordination starts from line 9

def interpret_car(lines):
    index = 5 
    ele_name  = lines[index].strip().split()
    index = 6 
    ele_num = [int(i) for i in lines[index].strip().split()]
    dict_car =  {ele_name[i]:ele_num[i] for i in range(0, len(ele_name))} 
    
    dict_car2 = {}
    for element in ele_name: 
        indice  = ele_name.index(element)
        n_start = sum(ele_num[0: indice+1]) - dict_car.get(element) + 1
        n_end   = sum(ele_num[0: indice+1]) 
        n_range = '%s---%s' %(n_start, n_end)
        dict_car2.update({element : n_range})
       
    dict_car3 = {y:x for x,y in dict_car2.iteritems()}
    return lines, ele_name, ele_num, dict_car, dict_car2, dict_car3

def determinelayers(z_cartesian):
    seq = sorted(z_cartesian)
    min = seq[0]
    layerscount = {}
    sets = [min]
    for j in range(len(seq)):
        if abs(seq[j]-min) >= options.thresthold:
            min = seq[j]
            sets.append(min)
    for i in range(1,len(sets)+1):
        layerscount[i] = []            
        if i > 1:
            for k in range(len(z_cartesian)):   
                if abs(z_cartesian[k]-sets[i-1]) <= options.thresthold and not abs(z_cartesian[k]-sets[i-2]) <= options.thresthold:
                    layerscount[i].append(k)            
        else:
            for k in range(len(z_cartesian)):   
                if abs(z_cartesian[k]-sets[i-1]) <= options.thresthold:
                    layerscount[i].append(k)
    return layerscount


def convert(start_num,**kwargs):
    tf = [] ;tf1=[]; tf2=[] 
    for i in range(start_num, num_atoms + start_num):
        x_cartesian.append(float(line[i].split()[0]) * a1 + float(line[i].split()[1]) * a2 + float(line[i].split()[2]) * a3)
        y_cartesian.append(float(line[i].split()[0]) * b1 + float(line[i].split()[1]) * b2 + float(line[i].split()[2]) * b3)
        z_cartesian.append(float(line[i].split()[0]) * z1 + float(line[i].split()[1]) * z2 + float(line[i].split()[2]) * z3)
         
        if len(line[i].split()) > 3:   # if  T T T exist, there are more than 3 elements in the list line[i].split()
            tf.append((line[i].split()[3]))
            tf1.append((line[i].split()[4]))
            tf2.append((line[i].split()[5]))
        else:
            tf.append(' ')   # if there is no T T T, use space instead. 
            tf1.append(' ')
            tf2.append(' ')
    #print layerscount
    if options.fixed and not  options.selected:
        layerscount =determinelayers(z_cartesian)
        fixedlayer = int(input("    Found %d layers, choose how many layers to be fixed------>" %len(layerscount)))
        #tf1=copy.deepcopy(tf);tf2=copy.deepcopy(tf)
        for i in range(1,len(layerscount)+1):
            if i <= fixedlayer: 
                for j in layerscount[i]:
                    tf[j] = " F " ; tf1[j] = " F " ;tf2[j] = " F "
            else:
                for k in layerscount[i]:
                    tf[k] = " T ";tf1[k] = " T " ;tf2[k] = " T "
    elif options.fixed and options.selected:
        #tf1=copy.deepcopy(tf);tf2=copy.deepcopy(tf)
        if kwargs:
            for i in kwargs['a'][0]:
                tf[i-1] = " F " ;tf1[i-1] = " F ";tf2[i-1] = " F "
    elif options.relax and options.selected:
        #tf1=copy.deepcopy(tf);tf2=copy.deepcopy(tf)
        TTT=''.join(kwargs['a'][1])
        if kwargs:
            for i in kwargs['a'][0]:
                tf[i-1] = " %s "  %TTT[0];tf1[i-1] = " %s " %TTT[1];tf2[i-1] = " %s " %TTT[2]
    else:
        tf1=copy.deepcopy(tf);tf2=copy.deepcopy(tf)

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
        %(x_cartesian[i], y_cartesian[i], z_cartesian[i], tf[i], tf1[i], tf2[i]))
    
    file_out.close()
    print '\n    -----------------------------------------------------------'
    print '    POSCAR with Cartesian Coordiations is named as %s_C' %(file_to_be_converted)
    if options.fixed:
        print '\n    Atom fix completed, if layers are not separated properly,\n    Please alter threstholds with flags -y or type -h for help!'
    print '    -----------------------------------------------------------\n'



def exe(*args):

    if line[7].strip().upper().startswith('S'):  # # With Selected T T T, coordination starts from line 9
        start_num = 9

        if  line[8].strip().upper().startswith('D'):
            print """
    This POSCAR has Direct Coordinations, Conversion is starting....
                  """
            if options.fixed and not  options.selected:
                print '\n    -----------------------------------------------------------------------\n'
                print u"    Now converting %s to cartesian coordinates\n    Then type how many layers to be fixed, from bottom to top." %options.files
                print '    -----------------------------------------------------------------------\n'
                convert(start_num)
            elif options.fixed and options.selected:
                print '\n    -----------------------------------------------------------------------\n'
                print u"    Now converting %s to cartesian coordinates\n    Then atoms selected will be fixed." %options.files
                print '    -----------------------------------------------------------------------\n'
                convert(start_num,a=args)
            elif options.relax and options.selected:
                print '\n    -----------------------------------------------------------------------\n'
                print u"    Now converting %s to cartesian coordinates\n    Then atoms selected will be relaxed." %options.files
                print '    -----------------------------------------------------------------------\n'                
                convert(start_num,a=args)
            else:
                convert(start_num)
            
        
        elif  line[8].strip().upper().startswith('C'):
            if  options.fixed and not  options.selected:    
                print '\n    -----------------------------------------------------------'
                print u"    Cartesian Coordinates found, only for fixing atoms!\n    Then type how many layers to be fixed, from bottom to top." 
                print '    -----------------------------------------------------------\n'
                convert(start_num)
            elif options.fixed and options.selected:
                print '\n    -----------------------------------------------------------------------\n'
                print u"    Cartesian Coordinates found, only for fixing atoms!\n    Then atoms selected will be fixed." 
                print '    -----------------------------------------------------------------------\n'                
                convert(start_num,a=args) 
            elif options.relax and options.selected:
                print '\n    -----------------------------------------------------------------------\n'
                print u"    Cartesian Coordinates found, only for relaxing atoms!\n    Then atoms selected will be relaxed." 
                print '    -----------------------------------------------------------------------\n'  
                convert(start_num,a=args)              
            else:
                print '\n    -----------------------------------------------------------'
                print "    This POSCAR has Cartesian Coordinations! Process is aborted!"
                print '    -----------------------------------------------------------\n'
                

    else : 
        print """
    -----------------------------------------------------------
    Pay Attetion! There is no TTT in coordinations part!
    -----------------------------------------------------------
    """
        
        start_num = 8 # without Selected, No  T T T , coordination starts from line 8 
        
        if  line[7].strip().upper().startswith('D'):
            print """
    This POSCAR has Direct Coordinations, Conversion starts....
    """
            if options.fixed and not options.selected:
                print '\n    -----------------------------------------------------------'
                print u"    Now converting %s to cartesian coordinates\n    Then type how many layers to be fixed, from bottom to top." %options.files
                print '    -----------------------------------------------------------\n'
            elif options.fixed and options.selected:
                print '\n    -----------------------------------------------------------'
                print u"    No TTT OR FFF,please add them beforehand\n    Now converting %s to cartesian coordinates\n    Then type how many layers to be fixed, from bottom to top." %options.files
                print '    -----------------------------------------------------------\n' 
                options.selected = 0 #deselect
            else:
                pass               
            convert(start_num)
        
        elif  line[7].strip().upper().startswith('C'):
            if  options.fixed and not options.selected:
                print '\n    -----------------------------------------------------------'
                print u"    Cartesian Coordinates found, only for fixing atoms!\n    Then type how many layers to be fixed, from bottom to top." 
                print '    -----------------------------------------------------------\n'
                convert(start_num)
            elif options.fixed and options.selected:
                print '\n    -----------------------------------------------------------'
                print u"    No TTT OR FFF,please add them beforehand\n    Cartesian Coordinates found, only for fixing atoms!\n    Then type how many layers to be fixed, from bottom to top." 
                print '    -----------------------------------------------------------\n'
                options.selected = 0 #deselect
                convert(start_num)             
            else:
                print '\n    -----------------------------------------------------------'
                print "    This POSCAR has Cartesian Coordinations! Process is aborted!"
                print '    -----------------------------------------------------------\n'

    file_read.close()

def get_atoms(selects,dict_car2,dict_car): # return a list of atoms which user selected
    ele_list = []
    for i in selects:
        if '-' not in i:   # for example: atoms 1-5 
            if i.isdigit():
                ele_list.append(int(i))
            elif i in dict_car2.keys(): # For Example: Ru C H O                                
                n_s, n_e = dict_car2.get(i).split('---')[:]
                ele_list.extend(range(int(n_s), int(n_e)+1)) 
            elif i == 'all':
                ele_list = range(1,sum([value for value in dict_car.values()])+1)
                break
            else:
                print('\n    Wrong Atom Selection!') 
        else:
            num_start = int(i.split('-')[0])
            num_end   = int(i.split('-')[1]) + 1
            ele_list.extend(range(num_start, num_end))
    return  ele_list 

def interpret_selected(s,selects=options.selected,arg_l = sys.argv[:]):
    sel=[]
    if selects:
        begin =  arg_l.index(s)
    for i in range(begin+1,len(arg_l)):
        if re.search(r'-[a-zA-Z]', arg_l[i]):
            break
        else:
            sel.append(arg_l[i]) 
    return sel

if __name__=='__main__':
    #lines, ele_name, ele_num, dict_car, dict_car2, dict_car3 = interpret_car(line)
    if options.relax:
        sels = interpret_selected('-r',options.relax,arg_l = sys.argv[:])   
        if not sels:
            sels=['TTT']
        print('\n    Your selective choice is ****  %s  ****!'%''.join(sels)) 
    if options.selected:
        lines, ele_name, ele_num, dict_car, dict_car2, dict_car3 = interpret_car(line)
        sel = interpret_selected('-s',options.selected,arg_l = sys.argv[:])
        ele_list = get_atoms(sel,dict_car2,dict_car)
        if not options.relax:
            exe(ele_list)
        else:
            exe(ele_list,sels)
    else:
        exe()
    
    #print dict_car, dict_car2, dict_car3
        