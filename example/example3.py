#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright 2011 by Paweł T. Jochym <pawel.jochym@ifj.edu.pl>
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
Example of elastic tensor calculation using VASP calculator and 
elastic module.
'''

from ase.lattice.spacegroup import crystal
from parcalc import ClusterVasp, ParCalculate
from elastic import Crystal
import ase.units as units
from numpy import array, linspace
import matplotlib.pyplot as plt
import numpy

# Create the MgO crystal
a = 4.194
cryst = Crystal(crystal(['Mg', 'O'], 
                [(0, 0, 0), (0.5, 0.5, 0.5)], 
                spacegroup=225,
                cellpar=[a, a, a, 90, 90, 90]))

# Create the calculator running on one, eight-core node.
# This is specific to the setup on my cluster.
# You have to adapt this part to your environment
calc=ClusterVasp(nodes=1,ppn=8)

# Assign the calculator to the crystal
cryst.set_calculator(calc)

# Set the calculation parameters
calc.set(prec = 'Accurate', 
            xc = 'PBE', 
            lreal = False, 
            nsw=30,
            ediff=1e-8, 
            ibrion=2,
            kpts=[3,3,3])

# Set the calculation mode first.
# Full structure optimization in this case.
# Not all calculators have this type of internal minimizer!
calc.set(isif=3)

print "Running initial optimization ... ",
print "Residual pressure: %.3f bar" % (
            cryst.get_isotropic_pressure(cryst.get_stress()))

# Clean up the directory
calc.clean()

# Switch to cell IDOF optimizer
calc.set(isif=2)

# Elastic tensor by internal routine
Cij, Bij=cryst.get_elastic_tensor(n=5,d=0.2)
print "Cij (GPa):", Cij/units.GPa


# Now let us do it (only c11 and c12) by hand 
# with 3rd order polynomial fitting the points.
sys=[]
for d in linspace(-0.2,0.2,10):
    sys.append(cryst.get_cart_deformed_cell(axis=0,size=d))
r=ParCalculate(sys,cryst.calc)
ss=[]
for s in r:
    ss.append([s.get_strain(cryst), s.get_stress()])

# Plot strain-stress relation
ss=[]
for p in r:
    ss.append([p.get_strain(cryst),p.get_stress()])
ss=array(ss)
lo=min(ss[:,0,0])
hi=max(ss[:,0,0])
mi=(lo+hi)/2
wi=(hi-lo)/2
xa=linspace(mi-1.1*wi,mi+1.1*wi, 50)
plt.plot(ss[:,0,0],ss[:,1,0],'k.')
plt.plot(ss[:,0,0],ss[:,1,1],'r.')

plt.axvline(0,ls='--')
plt.axhline(0,ls='--')


# C11 component
f=numpy.polyfit(ss[:,0,0],ss[:,1,0],3)
c11=f[-2]/units.GPa
plt.plot(xa,numpy.polyval(f,xa),'b-')

# C12 component
f=numpy.polyfit(ss[:,0,0],ss[:,1,1],3)
c12=f[-2]/units.GPa
plt.plot(xa,numpy.polyval(f,xa),'g-')

print 'C11 = %.3f GPa, C12 = %.3f GPa => K= %.3f GPa' % (c11, c12, (c11+2*c12)/3)

plt.show()

