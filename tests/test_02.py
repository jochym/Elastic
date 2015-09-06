# coding: utf-8
from __future__ import division, print_function

from ase import Atom, Atoms
from parcalc import ParCalculate, ClusterVasp
import elastic
from ase import io
from ase.lattice.spacegroup import crystal
from pyspglib import spglib
from time import sleep
import os, sys

def wait_for_results(systems):
    for n,r in enumerate(systems):
        print(n+1,end='')
        while not r.get_calculator().calc_finished() :
            print('.',end='')
            sys.stdout.flush()
            sleep(10)
        print(end=' ')
        sys.stdout.flush()
    print()


a=4 ; c=crystal('Si',[(0,0,0)],spacegroup=221,cellpar=[a,a,a,90,90,90])
spglib.get_spacegroup(c)

block=False

calc = ClusterVasp(block=block)

calc.set(prec = 'Accurate', xc = 'PBE', lreal = False,
            nsw=30, ediff=1e-8, ibrion=2, kpts=[3,3,3])

calc.set(isif=3)

c.set_calculator(calc)

l=ParCalculate(c,calc,cleanup=False,block=block)

s = l[0]

wait_for_results(l)

print(s.get_forces(), s.get_cell(), s.get_stress(), s.get_pressure())

calc.set(isif=2)
sl=s.get_BM_EOS(mode='gener')
res=ParCalculate(sl,calc,block=block,prefix='BMEOS_')

wait_for_results(res)

print(s.get_BM_EOS(data=res))

deformations=s.get_elastic_tensor(mode='deformations')
results=ParCalculate(deformations,calc,block=block,prefix='Cij_')

wait_for_results(results)

Cij,Bij=s.get_elastic_tensor(mode='restart',systems=results)

print(Cij)

