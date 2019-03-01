#!/bin/bash
for i in `ls -l |grep -v ^d|awk '{print $9}'` ;
 do sed -i "s/\t/    /g" ${i} ;
 sed -i "s/\r//g"  ${i} ;
 done
