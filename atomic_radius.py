# -*- coding:utf-8
import sys
dicts={'H': 53, 'He': 31, 'Li': 167, 'Be': 112, 'B': 87, 'C': 67, 'N': 56, 'O': 48, 'F': 42, 'Ne': 38, 'Na': 190, 'Mg': 145, 'Al': 118, 'Si': 111, 'P': 98, 'S': 88, 'Cl': 79, 'Ar': 71, 'K': 243, 'Ca': 194, 'Sc': 184, 'Ti': 176, 'V': 171, 'Cr': 166, 'Mn': 161, 'Fe': 156, 'Co': 152, 'Ni': 149, 'Cu': 145, 'Zn': 142, 'Ga': 136, 'Ge': 125, 'As': 114, 'Se': 103, 'Br': 94, 'Kr': 88, 'Rb': 265, 'Sr': 219, 'Y': 212, 'Zr': 206, 'Nb': 198, 'Mo': 190, 'Tc': 183, 'Ru': 178, 'Rh': 173, 'Pd': 169, 'Ag': 165, 'Cd': 161, 'In': 156, 'Sn': 145, 'Sb': 133, 'Te': 123, 'I': 115, 'Xe': 108, 'Cs': 298, 'Ba': 253, 'Hf': 208, 'Ta': 200, 'W': 193, 'Re': 188, 'Os': 185, 'Ir': 180, 'Pt': 177, 'Au': 174, 'Hg': 171, 'Tl': 156, 'Pb': 154, 'Bi': 143, 'Po': 135, 'At': 127, 'Rn': 120}

if len(sys.argv)<=2:
    print("usage: python atomic_radius.py H C")
    sys.exit(0)
a1=dicts.get(sys.argv[1],0)
a2=dicts.get(sys.argv[2],0)
#print(a1,a2)
if a1==0 or a2==0:
    raise IOError("Elements not found!")
print("Atomic distances is %f Å" %((a1+a2)/100.0))

