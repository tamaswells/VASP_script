#!/bin/bash

vaspkit='/public/home/apclab/nxu/vasp-exe/vaspkit.0.73/bin/vaspkit'
(echo '105';echo $1;echo) | ${vaspkit} 
(echo '402';echo '1';echo '1.5';echo $2) | ${vaspkit}
mv POSCAR POSCAR.bak
mv POSCAR_fix POSCAR
echo "第一个参数为cif文件名，第二个参数是从底部固定的层数！"
 