#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#    Copyright 1998-2011 by Paweł T. Jochym <pawel.jochym@ifj.edu.pl>
#
#    This file is part of Elastic.

#    Elastic is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    Elastic is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Elastic.  If not, see <http://www.gnu.org/licenses/>.

'''

.. _elastic-mod:

Elastic Module
^^^^^^^^^^^^^^

This module depends on :ref:`par-calc-mod` for parallelisation of 
independent calculations.

Elastic is a module for calculation of :math:`C_{ij}` components of elastic 
tensor from the strain-stress relation.
 
The strain components here are ordered in standard way which is different
to ordering in previous versions of the code.

The ordering is: :math:`u_{xx}, u_{yy}, u_{zz}, u_{yz}, u_{xz}, u_{xy}`.

The general ordering of :math:`C_{ij}` components is (except for triclinic symmetry and taking into account customary names of constants - e.g. 
:math:`C_{16} \\rightarrow C_{14}`):

.. math::
   C_{11}, C_{22}, C_{33}, C_{12}, C_{13}, C_{23}, 
   C_{44}, C_{55}, C_{66}, C_{16}, C_{26}, C_{36}, C_{45}

The functions outside of the Crystal class define the symmetry of the
:math:`C_{ij}` matrix. The matrix is N columns by 6 rows where the columns
corespond to independent elastic constants of the given crystal, while the rows
corespond to the canonical deformations of a crystal. The elements are the
second partial derivatives of the free energy formula for the crystal written
down as a quadratic form of the deformations with respect to elastic constant
and deformation. 

*Note:*
The elements for deformations :math:`u_{xy}, u_{xz}, u_{yz}`
have to be divided by 2 to properly match the usual definition 
of elastic constants.

See: [LL]_ L.D. Landau, E.M. Lifszyc, "Theory of elasticity"

There is some usefull summary also at: 
`ScienceWorld <http://scienceworld.wolfram.com/physics/Elasticity.html>`_

Class description
"""""""""""""""""
'''

from __future__ import print_function, division, absolute_import

import re
import sys
import string


import ase.io
from ase.atoms import Atoms

try :
    # Try new release of spglib
    import spglib as spg
except ImportError :
    # Old naming scheme
    from pyspglib import spglib as spg
    
from scipy.linalg import norm, lstsq
from scipy import optimize
from numpy.linalg import inv
from numpy import dot, diag, ones, reshape, linspace, array, mean
from math import acos, pi, cos, sin, tan
import ase.units as units


def BMEOS(v,v0,b0,b0p):
    return (b0/b0p)*(pow(v0/v,b0p) - 1)

def ctg(x):
    return cos(x)/sin(x)

def csc(x):
    return 1/sin(x)


def regular(u):
    '''
    Equation matrix generation for the regular (cubic) lattice.
    The order of constants is as follows:

    .. math::
       C_{11}, C_{12}, C_{44}
    '''
    uxx, uyy, uzz, uyz, uxz, uxy = u[0],u[1],u[2],u[3],u[4],u[5]
    return array([
    [uxx,   uyy + uzz,      0],
    [uyy,   uxx + uzz,      0],
    [uzz,   uxx + uyy,      0],
    [0,             0,              2*uyz],
    [0,             0,              2*uxz],
    [0,             0,              2*uxy]])

def tetragonal(u):
    '''
    Equation matrix generation for the tetragonal lattice.
    The order of constants is as follows:

    .. math::
       C_{11}, C_{33}, C_{12}, C_{13}, C_{44}, C_{14}
    '''
    uxx, uyy, uzz, uyz, uxz, uxy = u[0],u[1],u[2],u[3],u[4],u[5]
    return array(
    [[uxx,   0,    uyy,  uzz,      0,      0],
     [uyy,   0,    uxx,  uzz,      0,      0],
     [0,     uzz,  0,    uxx+uyy,  0,      0],
     [0,     0,    0,    0,        0,      2*uxy],
     [0,     0,    0,    0,        2*uxz,  0],
     [0,     0,    0,    0,        2*uyz,  0]])
 

def orthorombic(u):
    '''
    Equation matrix generation for the orthorombic lattice.
    The order of constants is as follows:

    .. math::
       C_{11}, C_{22}, C_{33}, C_{12}, C_{13}, C_{23}, 
       C_{44}, C_{55}, C_{66}
    '''
    uxx, uyy, uzz, uyz, uxz, uxy = u[0],u[1],u[2],u[3],u[4],u[5]
    return array(
    [[uxx,    0,    0,  uyy,  uzz,    0,    0,    0,    0],
    [0,     uyy,    0,  uxx,    0,  uzz,    0,    0,    0],
    [0,       0,  uzz,    0,  uxx,  uyy,    0,    0,    0],
    [0,       0,    0,    0,    0,    0,2*uyz,    0,    0],
    [0,       0,    0,    0,    0,    0,    0,2*uxz,    0],
    [0,       0,    0,    0,    0,    0,    0,    0,2*uxy]])
 

def trigonal(u):
    '''
    The matrix is constructed based on the approach from L&L
    using auxiliary coordinates: :math:`\\xi=x+iy`, :math:`\\eta=x-iy`.
    The components are calculated from free energy using formula 
    introduced in :ref:`symmetry` with appropriate coordinate changes.
    The order of constants is as follows:

    .. math::
       C_{11}, C_{33}, C_{12}, C_{13}, C_{44}, C_{14}
    '''
    #TODO: Not tested yet. 
    #TODO: There is still some doubt about the :math:`C_{14}` constant. 
    uxx, uyy, uzz, uyz, uxz, uxy = u[0],u[1],u[2],u[3],u[4],u[5]
    return array(
    [[   uxx,   0,    uyy,     uzz,     0,   2*uxz],
     [   uyy,   0,    uxx,     uzz,     0,  -2*uxz],
     [     0, uzz,      0, uxx+uyy,     0,   0],
     [     0,   0,      0,       0, 2*uyz,  -4*uxy],
     [     0,   0,      0,       0, 2*uxz,   2*(uxx-uyy)],
     [ 2*uxy,   0, -2*uxy,       0,     0,  -4*uyz]])

def hexagonal(u):
    '''
    The matrix is constructed based on the approach from L&L
    using auxiliary coordinates: :math:`\\xi=x+iy`, :math:`\\eta=x-iy`.
    The components are calculated from free energy using formula 
    introduced in :ref:`symmetry` with appropriate coordinate changes.
    The order of constants is as follows:
    
    .. math::
       C_{11}, C_{33}, C_{12}, C_{13}, C_{44}
    '''
    #TODO: Still needs good verification
    uxx, uyy, uzz, uyz, uxz, uxy  = u[0],u[1],u[2],u[3],u[4],u[5]
    return array(
    [[   uxx,   0,    uyy,     uzz,     0   ],
     [   uyy,   0,    uxx,     uzz,     0   ],
     [     0, uzz,      0, uxx+uyy,     0   ],
     [     0,   0,      0,       0, 2*uyz   ],
     [     0,   0,      0,       0, 2*uxz   ],
     [ 2*uxy,   0, -2*uxy,       0,     0   ]])

def monoclinic(u):
    '''
    Monoclinic group, the ordering of constants is:
    
    .. math::
       C_{11}, C_{22}, C_{33}, C_{12}, C_{13}, C_{23}, 
       C_{44}, C_{55}, C_{66}, C_{16}, C_{26}, C_{36}, C_{45}
    '''
    
    uxx, uyy, uzz, uyz, uxz, uxy = u[0],u[1],u[2],u[3],u[4],u[5]
    return array(
    [[uxx,  0,  0,uyy,uzz,  0,    0,    0,    0,uxy,  0,  0,  0],
     [  0,uyy,  0,uxx,  0,uzz,    0,    0,    0,  0,uxy,  0,  0],
     [  0,  0,uzz,  0,uxx,uyy,    0,    0,    0,  0,  0,uxy,  0],
     [  0,  0,  0,  0,  0,  0,2*uyz,    0,    0,  0,  0,  0,uxz],
     [  0,  0,  0,  0,  0,  0,    0,2*uxz,    0,  0,  0,  0,uyz],
     [  0,  0,  0,  0,  0,  0,    0,    0,2*uxy,uxx,uyy,uzz,  0]])


def triclinic(u):
    '''
    Triclinic crystals. 
    
    *Note*: This was never tested on the real case. Beware!
    
    The ordering of constants is:
    
    .. math::
       C_{11}, C_{22}, C_{33},
       C_{12}, C_{13}, C_{23},
       C_{44}, C_{55}, C_{66},
       C_{16}, C_{26}, C_{36}, C_{46}, C_{56},
       C_{14}, C_{15}, C_{25}, C_{45}
    '''
    # Based on the monoclinic matrix and not tested on real case.
    # If you have test cases for this symmetry send them to the author.
    uxx, uyy, uzz, uyz, uxz, uxy = u[0],u[1],u[2],u[3],u[4],u[5]
    return array(
    [[uxx,  0,  0,uyy,uzz,  0,    0,    0,    0,uxy,  0,  0,  0,  0,uyz,uxz,  0,  0],
     [  0,uyy,  0,uxx,  0,uzz,    0,    0,    0,  0,uxy,  0,  0,  0,  0,  0,uxz,  0],
     [  0,  0,uzz,  0,uxx,uyy,    0,    0,    0,  0,  0,uxy,  0,  0,  0,  0,  0,  0],
     [  0,  0,  0,  0,  0,  0,2*uyz,    0,    0,  0,  0,  0,uxy,  0,uxx,  0,  0,uxz],
     [  0,  0,  0,  0,  0,  0,    0,2*uxz,    0,  0,  0,  0,  0,uxy,  0,uxx,uyy,uyz],
     [  0,  0,  0,  0,  0,  0,    0,    0,2*uxy,uxx,uyy,uzz,uyz,uxz,  0,  0,  0,  0]])
    


#def ParCalculate(systems,calc):
#    for s in systems:
#        s.set_calculator(calc.copy())
#    calc.ParallelCalculate(systems,properties=['stress'])
#    return systems

from parcalc import ParCalculate
    

class CrystalInitError(Exception):
    def __str__(self):
        return '''The Crystal class should NEVER be created by itself - it is intended as 
        a mix-in base class for the Atoms class. Thus this constructor 
        just prints the error message and bails out.'''

class Crystal(Atoms):
    '''Backward compatibility class. To be removed later.'''
    pass

class ElasticCrystal:
    '''
    Mixin extension of standard ASE Atoms class designed to handle specifics of the
    crystalline materials. This code should, in principle, be folded into the 
    Atoms class in the future. At this moment it is too early to think about it.
    Additionally there are some aspects of this code which may be difficult to
    harmonize with the principles of the Atoms class. I am sure it is better,
    for now to leave this as a separate extension class.

    Basically, this class provides set of functions concerned with derivation 
    of elastic properties using "finite deformation approach" 
    (see the documentation for physics background information).
    '''

    def __init__(self):
        '''
        Dummy constructor for the Crystal class.
        The class should NEVER be created by itself - it is intended as 
        a mix-in base class for the Atoms class. Thus this constructor 
        just prints the error message and bails out.
        '''
        print('''Crystal class is not intended to be used directly! 
            You should never call it constructor. Read the docs or just 
            import elastic module and enjoy the new functionality of 
            the Atoms class!. Since this program is not going to work
            anyway I am bailing out right now.
            ''')
        raise CrystalInitError

    def get_lattice_type(self):
        '''
        Find the symmetry of the crystal using spglib symmetry finder.
        Assign to sg_name i sg_nr members name of the space group and
        its number extracted from the result. Based on the group number
        identify also the lattice type (assigned to sg_type member) and
        the Bravais lattice of the crystal (assigned to bravais member).
        The returned value is the lattice type number.
        The lattice type numbers are 
        (see also Crystal.ls, the numbering starts from 1): 
        
        Triclinic (1), Monoclinic (2), Orthorombic (3), Tetragonal (4)
        Trigonal (5), Hexagonal (6), Cubic (7)
        '''
        # Table of lattice types and correcponding group numbers dividing
        # the ranges. See get_lattice_type method for precise definition.

        lattice_types=[
                [3,   "Triclinic"],
                [16,  "Monoclinic"],
                [75,  "Orthorombic"],
                [143, "Tetragonal"],
                [168, "Trigonal"],
                [195, "Hexagonal"],
                [231, "Cubic"]
            ]

        sg=spg.get_spacegroup(self)
        m=re.match('([A-Z].*\\b)\s*\(([0-9]*)\)',sg)
        self.sg_name=m.group(1)
        self.sg_nr=int(m.group(2))
        
        for n,l in enumerate(lattice_types) :
            if self.sg_nr < l[0] :
                lattice=l[1]
                lattype=n+1
                break
        self.sg_type=lattype
        self.bravais=lattice
        return lattype
        
    def get_bulk_modulus(self,n=5, lo=0.98, hi=1.02, recalc=False):
        '''
        Calculate bulk modulus using the Birch-Murnaghan equation of state
        data calculated by get_BM_EOS routine (see). 
        The returned bulk modulus is a :math:`B_0` coefficient of the B-M EOS.
        The arguments are the same as in BM EOS function.
        '''
        if self._calc is None:
            raise RuntimeError('Crystal object has no calculator.')

        if recalc or getattr(self,'bm_eos',None) is None :
            self.get_BM_EOS(n,lo,hi,recalc)
        self.bulk_modulus=self.bm_eos[1]
        return self.bulk_modulus
        
    def get_pressure(self,s=None):
        '''
        Return *external* isotropic (hydrostatic) pressure in ASE units.
        If the pressure is positive the system is under external pressure.
        This is a convenience function.
        '''
        if s is None :
            s=self.get_stress()
        return -mean(s[:3])
        
    def get_BM_EOS(self,n=5, lo=0.98, hi=1.02, recalc=False, cleanup=True, mode='full', data=None):
        """
        Calculate Birch-Murnaghan Equation of State for the crystal:
        
        .. math::
           P(V)= \\frac{B_0}{B'_0}\\left[
           \\left({\\frac{V}{V_0}}\\right)^{-B'_0} - 1
           \\right]
        
        using n single-point structures ganerated from the 
        crystal (self) by the scan_volumes method between lo and hi 
        relative volumes. The BM EOS is fitted to the computed points by 
        least squares method. The returned value is a list of fitted 
        parameters: :math:`V_0, B_0, B_0'` if the fit succeded. 
        If the fitting fails the RuntimeError('Calculation failed') is reised.
        The data from the calculation and fit is stored in the bm_eos and pv
        members for future reference.
        
        *Note:* For now you have to set up the calculator to properly 
        optimize the structure without changing the volume at each point.
        There will be a way to specify basic types of the calculator 
        minimization at the later stage.
        """
        if self._calc is None:
            raise RuntimeError('Crystal object has no calculator.')

        if getattr(self,'bm_eos',None) is None or recalc :
            # NOTE: The calculator should properly minimize the energy
            # at each volume by optimizing the internal degrees of freedom
            # in the cell and/or cell shape without touching the volume.
            # TODO: Provide api for specifying IDOF and Full optimization 
            #       calculators. Maybe just calc_idof and calc_full members?
            if data is not None : # analyse results of previous calc
                res=data
            elif mode=='full' : # Make blocking calc of everything
                res=ParCalculate(self.scan_volumes(lo,hi,n),self.get_calculator(),cleanup=cleanup)
            elif mode=='gener' : # generate data for separate calc
                return self.scan_volumes(lo,hi,n)
            else :
                print('Error: Unrecognized mode and no data. Read the docs!')
                return
            
        #for r in res :
        #print(r.get_volume(), self.get_pressure(), r.get_cell())

            pvdat=array([[r.get_volume(),
                            self.get_pressure(r.get_stress()),
                            norm(r.get_cell()[:,0]),
                            norm(r.get_cell()[:,1]),
                            norm(r.get_cell()[:,2])] for r in res])
            #print(pvdat)

            # Fitting functions
            fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]
            errfunc = lambda p, x, y: fitfunc(p, x) - y

            # Estimate the initial guess assuming b0p=1
            # Limiting volumes
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
            #print(p0)
            p1, succ = optimize.leastsq(errfunc, p0[:], args=(pvdat[:,0],pvdat[:,1]))
            if not succ :
                raise RuntimeError('Calculation failed')
            self.bm_eos=p1
            self.pv=pvdat
        return self.bm_eos

    def get_elastic_tensor(self, n=5, d=2, mode='full', systems=None):
        '''
        Calculate elastic tensor of the crystal.
        It is assumed that the crystal is converged and optimized
        under intended pressure/stress.
        The geometry and stress at the call point is taken as
        the reference point. No additional optimization will be run.
        It is also assumed that the calculator is set to pure IDOF optimization.
        The size of used finite deformation is passed in d parameter as a
        percentage relative deformation. The n parameter defines number of
        deformed structures used in the calculation.
        '''
        # TODO: Provide API to enforce calculator selection
        
        # Deformation look-up table
        # Perhaps the number of deformations for trigonal 
        # system could be reduced to [0,3] but better safe then sorry
        deform={
            "Cubic": [[0,3], regular],
            "Hexagonal": [[0,2,3,5], hexagonal],
            "Trigonal": [[0,1,2,3,4,5], trigonal],
            "Tetragonal": [[0,2,3,5], tetragonal],
            "Orthorombic": [[0,1,2,3,4,5], orthorombic],
            "Monoclinic": [[0,1,2,3,4,5], monoclinic],
            "Triclinic": [[0,1,2,3,4,5], triclinic]
        }
        
        self.get_lattice_type()
        # Decide which deformations should be used
        axis, symm=deform[self.bravais]
        
        if mode!='restart':
            # Generate deformations if we are not in restart mode
            systems=[]
            for a in axis :
                if a <3 : # tetragonal deformation
                    for dx in linspace(-d,d,n):
                        systems.append(self.get_cart_deformed_cell(axis=a, size=dx))
                elif a<6 : # sheer deformation (skip the zero angle)
                    for dx in linspace(d/10.0,d,n):
                        systems.append(self.get_cart_deformed_cell(axis=a, size=dx))
        
        # Just generate deformations for manual calculation
        if mode=='deformations' :
            return systems
            
        if mode!='restart' :
            # Run the calculation if we are not restarting
            r=ParCalculate(systems,self.get_calculator())
        else :
            r=systems
        
        ul=[]
        sl=[]
        p=self.get_pressure()
        for g in r:
            ul.append(g.get_strain(self))
            # Remove the pressure from the stress tensor
            sl.append(g.get_stress()-array([p,p,p,0,0,0]))
        eqm=array(map(symm,ul))
        #print(eqm[0].shape, eqm.shape)
        eqm=reshape(eqm,(eqm.shape[0]*eqm.shape[1],eqm.shape[2]))
        #print(eqm)
        slm=reshape(array(sl),(-1,))
        #print(eqm.shape, slm.shape)
        #print(slm)
        Bij = lstsq(eqm,slm)
        #print(Bij[0] / units.GPa)
        # Calculate elastic constants from Birch coeff.
        # TODO: Check the sign of the pressure array in the B <=> C relation
        if (symm == orthorombic):
            Cij = Bij[0] - array([-p,-p,-p, p, p, p,-p,-p,-p])
        elif (symm == tetragonal):
            Cij = Bij[0] - array([-p,-p, p, p,-p,-p])
        elif (symm == regular):
            Cij = Bij[0] - array([-p, p,-p])
        elif (symm == trigonal):
            Cij = Bij[0] - array([-p,-p,p,p,-p,p])
        elif (symm == hexagonal):
            Cij = Bij[0] - array([-p,-p,p,p,-p])
        elif (symm == monoclinic):
            #TODO: verify this pressure array
            Cij = Bij[0] - array([-p,-p,-p, p, p, p,-p,-p,-p, p, p, p, p])
        elif (symm == triclinic):
            #TODO: verify this pressure array
            Cij = Bij[0] - array([-p,-p,-p, p, p, p,-p,-p,-p, p, p, p, p, p, p, p, p, p])
        return Cij, Bij

    def scan_pressures(self, lo, hi, n=5):
        '''
        Scan the pressure axis from lo to hi (inclusive) 
        using B-M EOS as the volume predictor.
        Pressure (lo, hi) in GPa
        '''
        # Inverse B-M EOS to get volumes from pressures
        # This will work only in limited pressure range p>-B/B'.
        # Warning! Relative, the V0 prefactor is removed.
        invbmeos = lambda b, bp, x: array([pow(b/(bp*xv+b),1/(3*bp)) for xv in x])
        
        eos=self.get_BM_EOS()
        
        # Limit negative pressures to 90% of the singularity value.
        # Beyond this B-M EOS is bound to be wrong anyway.
        lo=max(lo,-0.9*eos[1]/eos[2])
        
        scale=(eos[0]/self.get_volume())*invbmeos(eos[1], eos[2], 
                                                    linspace(lo,hi,num=n))
        #print(scale)
        uc=self.get_cell()
        sys=[Atoms(self) for s in scale]
        for n, s in enumerate(scale):
            sys[n].set_cell(s*uc,scale_atoms=True)
        
        return sys
        
        
    def scan_volumes(self, lo, hi, n):
        '''
        Provide set of crystals along volume axis from lo to hi (inclusive).
        No volume cell optimization is performed. Bounds are specified as 
        fractions (1.10 = 10% increase).
        '''
        scale=linspace(lo,hi,num=n)
        uc=self.get_cell()
        sys=[Atoms(self) for s in scale]
        for n, s in enumerate(scale):
            sys[n].set_cell(s*uc,scale_atoms=True)
        
        return sys

        
    def get_vecang_cell(self, uc=None):
        '''
        Compute A,B,C, alpha,beta,gamma cell params
        from the unit cell matrix (uc) or self.
        Angles in radians.
        '''
        if uc is None :
            uc=self.get_cell()
        ucv=[uc[i,:]/norm(uc[i,:]) for i in range(3)]
        uca=[acos(dot(ucv[(i+1)%3],ucv[(i+2)%3])) for i in range(3)]
        return [norm(uc[i,:]) for i in range(3)] + uca
        
    def get_deformed_cell(self, axis=0, size=1):
        '''
        Return the cell (with atoms) deformed along one
        cell parameter (0,1,2 = a,b,c ; 3,4,5 = alpha,beta,gamma) by
        size percent or size degrees (axis/angles).
        '''
        cryst=Crystal(self)
        uc=self.get_cell()
        if axis < 3 :
            uc[axis,:]=(1+size/100.0)*uc[axis,:]
        else :
            (a,b,c,alp,bet,gam)=cryst.get_vecang_cell()
            d=array([0.0, 0.0, 0.0])
            d[axis-3]=pi*size/180
            (alp,bet,gam)=array((alp,bet,gam))+d
            t=1 - (ctg(bet)*ctg(gam)-cos(alp)*csc(bet)*csc(gam))**2;
            if t<0.0 :
                print('''
                The parameters (alpha,beta,gamma)=(%f,%f,%f) are probably 
                incorrect and lead to imaginary coordinates. 
                This range of parameters is unsupported by this program 
                (and is, let me say, very strange for a crystal).
                Cennot continue, bye.''' % (alp,bet,gam))
                raise ValueError
            else :
                uc=[[a,0.0,0.0],
                    [b*cos(gam), b*sin(gam), 0],
                    [c*cos(bet), 
                        c*(cos(alp)/sin(gam) - cos(bet)*ctg(gam)),
                        c*sin(bet)*sqrt(t)]]
        cryst.set_cell(uc, scale_atoms=True)
        #print(cryst.get_cell())
        #print(uc)
        return cryst

    def get_cart_deformed_cell(self, axis=0, size=1):
        '''
        Return the cell (with atoms) deformed along one 
        of the cartesian directions 
        (0,1,2 = x,y,z ; sheers: 3,4,5 = yz, xz, xy) by
        size percent.
        '''
        cryst=Crystal(self)
        uc=cryst.get_cell()
        l=size/100.0
        L=diag(ones(3))
        if axis < 3 :
            L[axis,axis]+=l
        else :
            if axis==3 :
                L[1,2]+=l
            elif axis==4 :
                L[0,2]+=l
            else :
                L[0,1]+=l
        uc=dot(uc,L)
        cryst.set_cell(uc, scale_atoms=True)
        #print(cryst.get_cell())
        #print(uc)
        return cryst
        
    def get_strain(self,refcell=None):
        '''
        Return the strain tensor in the Voight notation as a conventional 
        6-vector. The calculation is done with respect to the crystal 
        geometry passed in refcell parameter.
        '''
        if refcell is None :
            refcell=self
        du=self.get_cell()-refcell.get_cell()
        m=refcell.get_cell()
        m=inv(m)
        u=dot(m,du)
        u=(u+u.T)/2
        return array([u[0,0], u[1,1], u[2,2], u[2,1], u[2,0], u[1,0]])




