# -*- coding:utf-8 -*-

from __future__ import division
import numpy as np
from glob import glob
import sys
import os

desc = '''
author: nxu
version:1.0
前提：安装了python2/3和numpy库。
功能：分割P4VASP 生成的data文件，并积分DOS，计算出 band center 和 电子数。
使用方法：python analysis_dos.py your_data_file 或者 python analysis_dos.py  前提是你的datafile名字叫dat。
tamas@zju.edu.cn
'''

if len(sys.argv) <= 1:
    print(u"请在后面带上dat文件名，否则默认查找dat")
    if not os.path.isfile("dat"):
        print(u"Error,dat文件不存在！")
        exit()
    else:
        datafile = "dat"
else:
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

def integrate(data,index):
    ener_dos = data[:, 0] * (data[:, index])
    fenmu = np.trapz(data[:, index], data[:, 0])
    fenzi = np.trapz(ener_dos, data[:, 0])
    center = fenzi / fenmu
    dianzimidu = data[np.where(data[:, 0] <= 0)]
    dianzishu = np.trapz(dianzimidu[:, index], dianzimidu[:, 0])
    return center,dianzishu

for i in glob("*_.dat"):
    print ("_____________")
    print (u"正在处理%s" % i)
    data = np.loadtxt(i, dtype=np.float64)
    columns = data.shape[1]
    if columns <=1:
        print (u"你的数据只有一列！")
        exit()
    elif columns == 2:
        center,dianzishu = integrate(data,1)
        print(u"band center为：%10.6f eV"% center)
        print (u"电子数目为：%d" % int(dianzishu))
        print ("_____________")
    elif columns == 3:
        print (u"发现SPIN1和SPIN2.")
        center,dianzishu = integrate(data,1)
        print(u"SPIN1 的 band center为：%10.6f eV"% center)
        print (u"SPIN1 的电子数目为：%d" % int(dianzishu))
        center,dianzishu = integrate(data,2)
        print(u"SPIN2 的 band center为：%10.6f eV"% center)
        print (u"SPIN2 的电子数目为：%d" % int(dianzishu))
        print ("_____________")
    else:
        print(u"超过4列！")
        exit()
