#!/bin/bash
# author: nxu
# version:1.0
# To extract fermi and orbital information from OUTCAR
# To use it : bash extract-fermi.sh
# tamas@zju.edu.cn

python -c 'import re;print re.findall(r"(E-fermi.*?)------------------------",open("OUTCAR").read(),re.S)[-1].rstrip() '
