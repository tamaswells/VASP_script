#!/usr/bin/env python
# -*- coding: utf-8 -*-
# http://staff.ustc.edu.cn/~zqj/posts/NVT-MD/
# AUTHOR: Qijing Zheng,UTST
import argparse
import ase
import sys
import numpy as np
from ase.io import read, write


def nose_mass(temperature, ndof, t0, L):
    '''
    Suggested Q:
        Shuichi Nosé, J. Chem. Phys., 81, 511(1984). 
    input:
    temperaute: in unit of Kelvin
    ndof: No. of degrees of freedom
    t0: The oscillation time in fs
    L: the length of the first basis vector
    '''

    # Q in Energy * Times**2
    qtmp = (t0 * 1E-15 / np.pi / 2)**2 * \
        2 * ndof * ase.units.kB * temperature \
        * ase.units._e

    # Q in AMU * Angstrom**2
    Q = qtmp / ase.units._amu / (L * 1E-10)**2

    return Q


def cnt_dof(atoms):
    '''
    Count No. of Degrees of Freedom
    '''
    if atoms.constraints:
        from ase.constraints import FixAtoms, FixScaled, FixedPlane, FixedLine
        sflags = np.zeros((len(atoms), 3), dtype=bool)
        for constr in atoms.constraints:
            if isinstance(constr, FixScaled):
                sflags[constr.a] = constr.mask
            elif isinstance(constr, FixAtoms):
                sflags[constr.index] = [True, True, True]
            elif isinstance(constr, FixedPlane):
                mask = np.all(np.abs(np.cross(constr.dir, atoms.cell)) < 1e-5,
                              axis=1)
                if sum(mask) != 1:
                    raise RuntimeError(
                        'VASP requires that the direction of FixedPlane '
                        'constraints is parallel with one of the cell axis')
                sflags[constr.a] = mask
            elif isinstance(constr, FixedLine):
                mask = np.all(np.abs(np.cross(constr.dir, atoms.cell)) < 1e-5,
                              axis=1)
                if sum(mask) != 1:
                    raise RuntimeError(
                        'VASP requires that the direction of FixedLine '
                        'constraints is parallel with one of the cell axis')
                sflags[constr.a] = ~mask

        return np.sum(~sflags)
    else:
        return len(atoms) * 3 - 3


def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='Real POSCAR to calculate the No. of Degrees of freedom')
    arg.add_argument('-u', dest='unit', action='store', type=str,
                     default='fs',
                     choices=['cm-1', 'fs'],
                     help='Default unit of input frequency')
    arg.add_argument('-T', '--temperature', dest='temperature',
                     action='store', type=float,
                     default=300,
                     help='The temperature.')
    arg.add_argument('-f', '--frequency', dest='frequency',
                     action='store', type=float,
                     default=40,
                     help='The frequency of the temperature oscillation')

    return arg.parse_args(cml)


if __name__ == '__main__':
    arg = parse_cml_args(sys.argv[1:])

    if arg.unit == 'cm-1':
        THzToCm = 33.3564095198152
        t0 = 1000 * THzToCm / arg.frequency
    else:
        t0 = arg.frequency

    pos  = read(arg.poscar)
    L    = np.linalg.norm(pos.cell, axis=1)[0]
    ndof = cnt_dof(pos)
    Q    = nose_mass(arg.temperature, ndof, t0, L)

    print("SMASS = {}".format(Q))
    