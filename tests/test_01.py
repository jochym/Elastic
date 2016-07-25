# Test case for the code. Calculate few well-known crystals

from __future__ import print_function, division

import os
import sys
from numpy import linspace, array, arange
import numpy
from math import pow

from matplotlib.pyplot import plot, show, figure, draw, axvline, axhline
from pylab import *

from ase.lattice.spacegroup import crystal
from ase.visualize import view
from ase.units import GPa

import elastic
from parcalc import ParCalculate, ClusterVasp
from elastic import BMEOS

def banner(msg):
    print()
    print(60*'=')
    print(max(0,(60-len(msg))//2-1)*' ',msg)
    print(60*'=')

def secban(msg):
    print()
    print(max(0,(60-len(msg))//2-1)*' ',msg)
    print(60*'-')


banner('Structure optimization on MgO')

a = 4.291
MgO = crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[a, a, a, 90, 90, 90])
               
##################################
# Provide your own calculator here
##################################
calc=ClusterVasp(nodes=1,ppn=8)
# The calculator must be runnable in an isolated directory
# Without disturbing other running instances of the same calculator
# They are run in separate processes (not threads!)

MgO.set_calculator(calc)
calc.set(prec = 'Accurate', xc = 'PBE', lreal = False, isif=2, nsw=20, ibrion=2, kpts=[1,1,1])

print("Residual pressure: %.3f GPa" % (MgO.get_pressure()/GPa))
calc.clean()

banner('Volume scan on MgO (+/- 5%)')

systems=[]
for av in numpy.linspace(a*0.95,a*1.05,5):
    systems.append(crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
               cellpar=[av, av, av, 90, 90, 90]))

pcalc=ClusterVasp(nodes=1,ppn=8)
pcalc.set(prec = 'Accurate', xc = 'PBE', lreal = False, isif=2, nsw=20, ibrion=2, kpts=[1,1,1])
res=ParCalculate(systems,pcalc)

v=[]
p=[]
for s in res :
    v.append(s.get_volume())
    p.append(s.get_pressure()/GPa)

plot(v,p,'o')
show()
    


# You can specify the directory with prepared VASP crystal for the test run
# or run through all prepared cases.
if len(sys.argv)>1 :
    crystals=[crystal(ase.io.read(sys.argv[1]+'/CONTCAR'))]
else :
    # Pre-cooked test cases
    crystals=[]
    
    # Cubic
    a = 4.194
    crystals.append(crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)],
        spacegroup=225, cellpar=[a, a, a, 90, 90, 90]))
#        a = 4.194
#        crystals.append(Crystal(crystal(['Ti', 'C'], [(0, 0, 0), (0.5, 0.5, 0.5)],
#            spacegroup=225, cellpar=[a, a, a, 90, 90, 90])))
    # Tetragonal
    a = 4.60
    c = 2.96
    crystals.append(crystal(['Ti', 'O'], [(0, 0, 0), (0.302, 0.302, 0)],
        spacegroup=136, cellpar=[a, a, c, 90, 90, 90]))
    # Trigonal (this is metal - for sure the k spacing will be too small)
    a = 4.48
    c = 11.04
    crystals.append(crystal(['Sb'], [(0, 0, 0.24098)],
        spacegroup=166, cellpar=[a, a, c, 90, 90, 120]))


banner("Running tests for Elastic")
# Iterate over all crystals. 
# We do not paralelize over test cases for clarity.
for cryst in crystals[:] :

    cryst.get_lattice_type()
    banner('Calculating: %s, %s, %s' % 
            (cryst.get_chemical_formula(), cryst.bravais, cryst.sg_name))
    # Setup the calculator
    calc=ClusterVasp(nodes=1,ppn=8)
    cryst.set_calculator(calc)
    calc.set(prec = 'Accurate', 
                xc = 'PBE', 
                lreal = False, 
                isif=4, 
                nsw=30,
                ediff=1e-8, 
                ibrion=2)
    calc.set(kpts=[3,3,3])
    
    # Run the calculations
    
    # Optimize the structure first
    calc.set(isif=3)
    
    secban('Structure optimization')
    
    # Run the internal optimizer
    print("Residual pressure: %.3f GPa" % (
                cryst.get_pressure()/GPa))
    print("Residual stress (GPa):", cryst.get_stress()/GPa)

    calc.clean()
    cryst.get_lattice_type()

    secban('Lattice parameters')

    print(cryst.get_vecang_cell())
    print(cryst.bravais, cryst.sg_type, cryst.sg_name, cryst.sg_nr)
    
    #view(cryst)

    # Switch to cell shape+IDOF optimizer
    calc.set(isif=4)

    secban('Birch-Murnaghan EOS')

    # Calculate few volumes and fit B-M EOS to the result
    # Use +/-3% volume deformation and 5 data points
    fit=cryst.get_BM_EOS(n=5,lo=0.97,hi=1.03)
    
    # Get the P(V) data points just calculated
    pv=array(cryst.pv)
    
    # Sort data on the first column (V)
    pv=pv[pv[:,0].argsort()]
    
    # Print the fitted parameters
    print("V0=%.3f A^3 ; B0=%.2f GPa ; B0'=%.3f ; a0=%.5f A" % ( 
            fit[0], fit[1]/GPa, fit[2], pow(fit[0],1./3)))
            
    v0=fit[0]

    # B-M EOS for plotting
    fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]

    # Ranges - the ordering in pv is not guarateed at all!
    # In fact it may be purely random.
    x=array([min(pv[:,0]),max(pv[:,0])])
    y=array([min(pv[:,1]),max(pv[:,1])])

    
    # Plot b/a and c/a as a function of v/v0
    figure(1)
    plot(pv[:,0]/v0,pv[:,3]/pv[:,2],'o')
    plot(pv[:,0]/v0,pv[:,4]/pv[:,2],'x-')
    #print(pv[:,4]/pv[:,2])
    axvline(1,ls='--')
    title('b/a and c/a as a function of v/v0')
    draw()
    
    # Plot the P(V) curves and points for the crystal
    figure(2)
    # Plot the points
    plot(pv[:,0]/v0,pv[:,1],'o')
    
    # Mark the center P=0 V=V0
    axvline(1,ls='--')
    axhline(0,ls='--')

    # Plot the fitted B-M EOS through the points
    xa=linspace(x[0],x[-1],20)
    #xa=v0*linspace(0.90,1.10,20)
    plot(xa/v0,fitfunc(fit,xa),'-')
    title('Birch-Murnaghan Equation of state')
    draw()
    
    # Scan over deformations
    
    # Switch to IDOF optimizer
    calc.set(isif=2)

    secban('Elastic tensor by internal routine')

    # Elastic tensor by internal routine
    Cij, Bij=cryst.get_elastic_tensor(n=5,d=0.33)
    print("Cij (GPa):", Cij/GPa)
    
    calc.clean()
    
    secban('C11 and C12 by hand')
    
    # Now let us do it (only c11 and c12) by hand
    sys=[]
    for d in linspace(-0.5,0.5,6):
        sys.append(cryst.get_cart_deformed_cell(axis=0,size=d))
    r=ParCalculate(sys,cryst.calc)
    ss=[]
    for s in r:
        ss.append([s.get_strain(cryst), s.get_stress()])
    # Plot strain-stress relation
    figure(3)

    ss=[]
    for p in r:
        ss.append([p.get_strain(cryst),p.get_stress()])
    ss=array(ss)
    lo=min(ss[:,0,0])
    hi=max(ss[:,0,0])
    mi=(lo+hi)/2
    wi=(hi-lo)/2
    xa=linspace(mi-1.1*wi,mi+1.1*wi, 50)
    plot(ss[:,0,0],ss[:,1,0],'k.')
    plot(ss[:,0,0],ss[:,1,1],'r.')
    
    # C11 component
    f=numpy.polyfit(ss[:,0,0],ss[:,1,0],3)
    c11=f[-2]/GPa
    plot(xa,numpy.polyval(f,xa),'b-')
    
    # C12 component
    f=numpy.polyfit(ss[:,0,0],ss[:,1,1],3)
    c12=f[-2]/GPa
    plot(xa,numpy.polyval(f,xa),'g-')
    print('C11 = %.3f GPa, C12 = %.3f GPa => K= %.3f GPa (cubic only)' % (c11, c12, (c11+2*c12)/3))
    axvline(0,ls='--')
    axhline(0,ls='--')
    title('C11, C12 manual fits')
    draw()

#        # Now make a short scan over pressures
#        
#        # Switch just shape+IDOF optimization
#        calc.set(isif=4)
#        sys=[]
#        sys=cryst.scan_pressures(-5.0, 5.0, 3)
#        r=ParCalculate(sys,cryst.calc)
#        print("Pressure scan (GPa):",end=" ")
#        for s in r :
#            print(cryst.get_pressure(s.get_stress())/GPa, end=" ")
#        print()
#        vl=array([s.get_volume() for s in r])
#        pl=array([cryst.get_pressure(s.get_stress())/GPa for s in r])
#        figure(2)
#        plot(vl/v0,pl,'+')
#        
#        # Check for proper inverse eos
#        invbmeos = lambda b, bp, x: array([pow(b/(bp*xv+b),1/bp) for xv in x])
#        xa=linspace(max(-8,-0.9*fit[1]/fit[2]),8,20)
#        ya=invbmeos(fit[1],fit[2],xa)
#        plot(ya,xa,'-')

banner('Tests finished, showing plots')

show()



