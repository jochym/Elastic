#!/usr/bin/python

import re
import sys
import string


import ase.io
from ase.atoms import Atoms
from pyspglib import spglib as spg
from pvasp import ParCalculate

def BMEOS(v,v0,b0,b0p):
    return (b0/b0p)*(pow(v0/v,b0p) - 1)

class Crystal(Atoms):

    def __init__(self, *args, **kwargs):
        Atoms.__init__(self, *args, **kwargs)
        self.recalc_bulk=True
        self.bulk_modulus=None
        self.bm_eos=None
        
    ls=[
        [3,   "Triclinic"],
        [16,  "Monoclinic"],
        [75,  "Orthorombic"],
        [143, "Tetragonal"],
        [168, "Trigonal"],
        [195, "Hexagonal"],
        [231, "Cubic"]
    ]

    def get_lattice_type(self):
        sg=spg.get_spacegroup(self)
        m=re.match('([A-Z].*\\b)\s*\(([0-9]*)\)',sg)
        self.sg_name=m.group(1)
        self.sg_nr=string.atoi(m.group(2))
        
        for n,l in enumerate(Crystal.ls) :
            if self.sg_nr < l[0] :
                lattice=l[1]
                lattype=n+1
                break
        self.sg_type=lattype
        self.bravais=lattice
        return lattype
        
    def get_bulk_modulus(self,n=5, lo=0.98, hi=1.02, recalc=False):
        """Calculate bulk modulus."""
        if self._calc is None:
            raise RuntimeError('Crystal object has no calculator.')

        if self.bm_eos == None :
            self.get_BM_EOS(n,lo,hi,recalc)
        self.bulk_modulus=self.bm_eos[1]
        return self.bulk_modulus
        
    def get_BM_EOS(self,n=5, lo=0.98, hi=1.02, recalc=False):
        """Calculate Birch-Murnaghan EOS."""
        if self._calc is None:
            raise RuntimeError('Crystal object has no calculator.')

        scale=linspace(lo,hi,num=n)
        if self.bulk_modulus==None or recalc :
            uc=self.get_cell()
            sys=[Atoms(self) for s in scale]
            for n, s in enumerate(scale):
                sys[n].set_cell(s*uc,scale_atoms=True)
            
            res=ParCalculate(sys,self.calc)
            
            # Volume in A^3 and pressure in GPa
            pvdat=array([[r.get_volume(),1e-4*r.get_isotropic_pressure(r.get_stress())] for r in res])
            #print pvdat

            # Fitting functions
            fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]
            errfunc = lambda p, x, y: fitfunc(p, x) - y

            # Estite the initial guess assuming b0p=1
            v1=min(pvdat[:,0])
            v2=max(pvdat[:,0])
            # The pressure is falling with the growing volume
            p2=min(pvdat[:,1])
            p1=max(pvdat[:,1])
            b0=(p1*v1-p2*v2)/(v2-v1)
            v0=v1*(p1+b0)/b0
            # Initial guess
            p0=[v0,b0,1]
            #Fitting
            #print p0
            p1, succ = optimize.leastsq(errfunc, p0[:], args=(pvdat[:,0],pvdat[:,1]))
            if not succ :
                raise RuntimeError('Calculation failed')
            self.bm_eos=p1
            self.pv=pvdat
        return self.bm_eos





if __name__ == '__main__':
    import os
    from numpy import linspace, array, arange
    from scipy import stats, optimize
    from math import pow

    from matplotlib.pyplot import plot,show
    from ase.lattice.spacegroup import crystal
    from pvasp import ClusterVasp

    if len(sys.argv)>1 :
        cryst=Crystal(ase.io.read(sys.argv[1]+'/CONTCAR'))
    else :
        a = 4.194
        MgO = crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
                cellpar=[a, a, a, 90, 90, 90])
        cryst=Crystal(MgO)

    cryst.get_lattice_type()

    print cryst.bravais, cryst.sg_type, cryst.sg_name, cryst.sg_nr

    calc=ClusterVasp(nodes=1,ppn=8)
    cryst.set_calculator(calc)
    calc.set(prec = 'Accurate', xc = 'PBE', lreal = False, isif=2, nsw=20, ibrion=2)

    cryst.calc.set(kpts=[5,5,5])
    fit=cryst.get_BM_EOS(n=10)
    pv=cryst.pv[:]
    print "V0=%.3f B0=%.2f B0'=%.3f a0=%.5f" % ( fit[0], fit[1], fit[2], pow(fit[0],1./3) )
    print "B0=", cryst.get_bulk_modulus()

    fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]

    plot(pv[:,0],pv[:,1],'o')
    x=array([pv[0,0],pv[-1,0]])
    y=array([pv[0,1],pv[-1,1]])
    plot([fit[0],fit[0]],y,'--')
    xa=linspace(x[0],x[-1],20)
    plot(xa,fitfunc(fit,xa),'-')
    show()




