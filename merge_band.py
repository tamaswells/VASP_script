# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
from glob import glob
import sys
import os

desc = '''
author: nxu
version:1.1
前提：安装了python2/3和numpy库。
功能：分割P4VASP 生成的data文件，合并多列y坐标，导入origin绘制能带图。
使用方法：python merge_band.py your_data_file 
tamas@zju.edu.cn
'''

if len(sys.argv) <= 1:
    print(u"请在后面带上data文件名，否则默认查找data文件")
    if not os.path.isfile("data"):
        print(u"Error,data文件不存在！")
        exit()
    else:
        print(u"正在处理data文件")
        datafile = "data"
        #flag = 1
else:
    if len(sys.argv) == 2:
        print(u"正在处理data文件")
        #flag = 1
    else:
        print(u"第3个参数以后将会舍去~")
        flag = 3
    datafile = sys.argv[1]

try:
    with open(datafile, 'r') as reader:
        datas = reader.read()
except:
    print(u"你指定的文件不存在！")
    exit()

index = 0
for i in datas.split("\n\n"):
    if len(i) < 100:
        continue
    filename = "dos_" + str(index) + "_.dat"
    with open(filename, 'w') as writer:
        writer.writelines(i)
    index += 1

index = 0 
for i in glob("*_.dat"):    
    datas = np.loadtxt(i, dtype=np.float64)    
    if index == 0:
        temp = datas
    else:
        #print datas[:,1].shape,temp.shape
        temp=np.c_[temp,datas[:,1]]
    index+=1
    columns = datas.shape[1]
#print temp
np.savetxt("new.dat",temp,fmt="%.4f",delimiter=" ")
for i in glob("*_.dat"):
    os.remove(i)
print (u"已完成")
