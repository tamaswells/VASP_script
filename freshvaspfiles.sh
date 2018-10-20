#!/bin/bash
shopt -s extglob #enable extglob (Extended Pattern Matching)
rm -rf !(INCAR|KPOINTS|POTCAR|POSCAR|vasp.pbs) #only keep these files