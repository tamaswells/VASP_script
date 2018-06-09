#!/bin/bash
#modified from VASP tutorial's script By nxu


if [ -z $1 ] ; then
    outcar="OUTCAR"
	echo 'Default OUTCAR'
else
    outcar=$1
fi

homo=`awk '/NELECT/ {print $3/2}' $outcar`
lumo=`awk '/NELECT/ {print $3/2+1}' $outcar`
nkpt=`awk '/NKPTS/ {print $4}' $outcar`

e1=`grep "     $homo     " $outcar | head -$nkpt | sort -n -k 2 | tail -1 | awk '{print $2}'`
e2=`grep "     $lumo     " $outcar | head -$nkpt | sort -n -k 2 | head -1 | awk '{print $2}'`

echo "HOMO: band:" $homo " E=" $e1
echo "LUMO: band:" $lumo " E=" $e2
bandgap=$(echo "scale=4;${e2}-${e1}" | bc)
echo "Bandgap: E=${bandgap}eV"
