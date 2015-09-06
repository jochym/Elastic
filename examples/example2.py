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
Example of Birch-Murnaghan EOS calculation using VASP calculator and 
elastic module.
'''

from ase.lattice.spacegroup import crystal
from parcalc import ClusterVasp, ParCalculate
from elastic import Crystal, BMEOS
import ase.units as units
from numpy import array, linspace
import matplotlib.pyplot as plt

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

# Lets extract optimized lattice constant.
# MgO is cubic so a is a first diagonal element of lattice matrix
a=cryst.get_cell()[0,0]

# Clean up the directory
calc.clean()

# Switch to cell shape+IDOF optimizer
calc.set(isif=4)

# Calculate few volumes and fit B-M EOS to the result
# Use +/-3% volume deformation and 5 data points
fit=cryst.get_BM_EOS(n=5,lo=0.97,hi=1.03)

# Get the P(V) data points just calculated
pv=array(cryst.pv)

# Sort data on the first column (V)
pv=pv[pv[:,0].argsort()]

# Print just fitted parameters
print "V0=%.3f A^3 ; B0=%.2f GPa ; B0'=%.3f ; a0=%.5f A" % ( 
        fit[0], fit[1]/units.GPa, fit[2], pow(fit[0],1./3))
        
v0=fit[0]

# B-M EOS for plotting
fitfunc = lambda p, x: [BMEOS(xv,p[0],p[1],p[2]) for xv in x]

# Ranges - the ordering in pv is not guarateed at all!
# In fact it may be purely random.
x=array([min(pv[:,0]),max(pv[:,0])])
y=array([min(pv[:,1]),max(pv[:,1])])


# Plot the P(V) curves and points for the crystal
# Plot the points
plt.plot(pv[:,0]/v0,pv[:,1],'o')

# Mark the center P=0 V=V0
plt.axvline(1,ls='--')
plt.axhline(0,ls='--')

# Plot the fitted B-M EOS through the points
xa=linspace(x[0],x[-1],20)
plt.plot(xa/v0,fitfunc(fit,xa),'-')
plt.show()



