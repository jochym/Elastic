# coding: utf-8
from __future__ import division, print_function

from ase import Atom, Atoms
import parcalc
import elastic
from ase import io
from ase.lattice.spacegroup import crystal
from pyspglib import spglib
from time import sleep
import os

a=4 ; c=crystal('Si',[(0,0,0)],spacegroup=221,cellpar=[a,a,a,90,90,90])
spglib.get_spacegroup(c)

calc = parcalc.ClusterVasp(block=False)

calc.set(prec = 'Accurate', xc = 'PBE', lreal = False,
            nsw=30, ediff=1e-8, ibrion=2, kpts=[3,3,3])

calc.set(isif=3)

c.set_calculator(calc)

print('Do the calc', os.getcwd())

l=parcalc.ParCalculate(c,calc,cleanup=False,block=True)

s = l[0]
while not s.get_calculator().calc_finished() :
    print('.',end='')
    sleep(10)

print(s.get_forces(), s.get_cell(), s.get_stress(), s.get_pressure())
