#!/usr/bin/python

import re
import sys
import string


import ase.io
from ase.atoms import Atoms
from pyspglib import spglib as spg
from pvasp import ParCalculate
from scipy.linalg import norm
from numpy import dot
from math import acos, pi

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
            pvdat=array([[r.get_volume(),
                            1e-4*r.get_isotropic_pressure(r.get_stress()),
                            norm(r.get_cell()[:,0]),
                            norm(r.get_cell()[:,1]),
                            norm(r.get_cell()[:,2])] for r in res])
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

    def get_elastic_tensor(self):
        pass
    
    def get_vecang_cell(self, uc=None):
        '''
        Compute A,B,C, alpha,beta,gamma cell params
        from the unit cell matrix (uc)
        '''
        if uc==None :
            uc=self.get_cell()
        ucv=[uc[i,:]/norm(uc[i,:]) for i in range(3)]
        uca=[180*acos(dot(ucv[(i+1)%3],ucv[(i+2)%3]))/pi for i in range(3)]
        return [norm(uc[i,:]) for i in range(3)] + uca
        
    def get_deformed_cell(self, axis=0, size=1):
        '''
        Return the cell (with atoms) deformed along one
        cell parameter (0,1,2 = a,b,c ; 3,4,5 = alpha,beta,gamma) by
        size percent or size degrees (axis/angles).
        '''
        cryst=Crystal(self)
        uc=cryst.get_cell()
        if axis < 3 :
            uc[axis,:]=(1+size/100.0)*uc[axis,:]
            cryst.set_cell(uc, scale_atoms=True)
            #print cryst.get_cell()
            #print uc
        else :
            pass
        return cryst




if __name__ == '__main__':
    import os
    from numpy import linspace, array, arange
    import numpy
    from scipy import stats, optimize
    from math import pow

    from matplotlib.pyplot import plot, show, figure, draw, axvline, axhline
    from ase.lattice.spacegroup import crystal
    from ase.visualize import view
    from pvasp import ClusterVasp

    if len(sys.argv)>1 :
        crystals=[Crystal(ase.io.read(sys.argv[1]+'/CONTCAR'))]
    else :
        crystals=[]
        
        # Cubic
        a = 4.194
        crystals.append(Crystal(crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)],
            spacegroup=225, cellpar=[a, a, a, 90, 90, 90])))
        # Tetragonal
        a = 4.59
        c = 2.96
        crystals.append(Crystal(crystal(['Ti', 'O'], [(0, 0, 0), (0.302, 0.302, 0)],
            spacegroup=136, cellpar=[a, a, c, 90, 90, 90])))
        # Trigonal
        a = 4.38
        c = 11.55
        crystals.append(Crystal(crystal(['Sb'], [(0, 0, 0.24098)],
            spacegroup=166, cellpar=[a, a, c, 90, 90, 120])))



    for cryst in crystals :

        cryst.get_lattice_type()

        print cryst.get_vecang_cell()
        print cryst.bravais, cryst.sg_type, cryst.sg_name, cryst.sg_nr
        
        view(cryst)
        calc=ClusterVasp(nodes=1,ppn=8)
        cryst.set_calculator(calc)
        calc.set(prec = 'Accurate', 
                    xc = 'PBE', 
                    lreal = False, 
                    isif=4, 
                    nsw=30, 
                    ibrion=2)

        cryst.calc.set(kpts=[3,3,3])
        fit=cryst.get_BM_EOS()
        pv=array(cryst.pv)
        
        # Sort data on the first column (V)
        pv=pv[pv[:,0].argsort()]
        
        print "V0=%.3f B0=%.2f B0'=%.3f a0=%.5f" % ( 
                fit[0], fit[1], fit[2], pow(fit[0],1./3) )
                
        v0=fit[0]

        fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]

        x=array([min(pv[:,0]),max(pv[:,0])])
        y=array([min(pv[:,1]),max(pv[:,1])])

        figure(1)
        plot(pv[:,0]/v0,pv[:,2],'x')
        plot(pv[:,0]/v0,pv[:,3],'+')
        plot(pv[:,0]/v0,pv[:,4],'o')
        axvline(1,ls='--')
        draw()
        
        figure(3)
        plot(pv[:,0]/v0,pv[:,3]/pv[:,2],'o')
        plot(pv[:,0]/v0,pv[:,4]/pv[:,2],'x-')
        #print pv[:,4]/pv[:,2]
        axvline(1,ls='--')
        draw()
        
        figure(2)
        plot(pv[:,0]/v0,pv[:,1],'o')
        axvline(1,ls='--')
        axhline(0,ls='--')
        xa=linspace(x[0],x[-1],20)
        plot(xa/v0,fitfunc(fit,xa),'-')
        draw()
        
        cryst.get_deformed_cell(axis=2,size=2)
    
    show()




