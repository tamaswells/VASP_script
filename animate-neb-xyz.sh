#!/bin/bash
for i in $(echo `ls -F | grep -E  "^[0-9]+\/$" | tr -d '/'` | sort) ; do cat ${i}/POSCAR.xyz >> total.xyz ; done