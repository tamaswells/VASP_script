#/bin/bash

# author: nxu
# version:1.0
# To determinr whether entropy T*S is less than 1mev
# To use it : bash sigma.sh
# tamas@zju.edu.cn

aveenergy=$(echo "scale=5;$(grep "entropy T\*S" OUTCAR | tail -n 1 |  awk '{print $5}')/$(sed -n 7p POSCAR | tr -d "\r"  |  awk '{ for(i=1;i<=NF;i++) sum+=$i; print sum}')" | bc)
sigma=`grep "SIGMA" INCAR | awk '{split($0,a,"=" ); print a[2]}'`
absaveenergy=${aveenergy#-}
if [ `echo "$absaveenergy < 0.001" | bc` -eq 1 ] ; then
	echo "Your sigma $sigma is okay; for average T*S $absaveenergy is < 0.001 eV "
else
	echo "Your sigma $sigma is BAD; for average T*S $absaveenergy  is > 0.001 eV,try to decrease sigma!"
fi	

