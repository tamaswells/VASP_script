#!/bin/bash

#!/usr/bin/env bash
# Create a GGA_PAW POTCAR file byconcatenation of POTCAR files
# BigBro 2017-07-31 TGN
# To Use it： potcar.sh Cu C H O
# Define local potpaw_GGA pseudopotentialrepository:
repo="/data1/potentials"
# Check if older version of POTCAR ispresent
if [ -f POTCAR ] ; then
	mv -f POTCAR old-POTCAR
echo " ** Warning: old POTCAR file found and renamed to 'old-POTCAR'."
fi
# Main loop - concatenate the appropriatePOTCARs (or archives)
for i in $*
do
if test -f $repo/$i/POTCAR ; then
	cat $repo/$i/POTCAR>> POTCAR
elif test -f $repo/$i/POTCAR.Z ; then
	zcat $repo/$i/POTCAR.Z >> POTCAR
elif test -f $repo/$i/POTCAR.gz ; then
	gunzip -c $repo/$i/POTCAR.gz >> POTCAR
else
	echo $i
	echo "Warning: No suitable POTCAR for element, Skipped this element."
fi
done
