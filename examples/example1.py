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
Example of parallel volume scan using VASP calculator and parcalc module.
'''

from ase.lattice.spacegroup import crystal
from parcalc import ClusterVasp, ParCalculate
import ase.units as units
import numpy
import matplotlib.pyplot as plt

a = 4.194
cryst = crystal(['Mg', 'O'], 
                [(0, 0, 0), (0.5, 0.5, 0.5)], 
                spacegroup=225,
                cellpar=[a, a, a, 90, 90, 90])

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

sys=[]
# Iterate over lattice constant in the +/-5% range
for av in numpy.linspace(a*0.95,a*1.05,5):
    sys.append(crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)], 
                spacegroup=225, cellpar=[av, av, av, 90, 90, 90]))
                   
# Define the template calculator for this run
# We can use the calc from above. It is only used as a template.
# Just change the params to fix the cell volume
calc.set(isif=2)

# Run the calculation for all systems in sys in parallel
# The result will be returned as list of systems res
res=ParCalculate(sys,calc)

# Collect the results
v=[]
p=[]
for s in res :
    v.append(s.get_volume())
    p.append(s.get_isotropic_pressure(s.get_stress()))

# Plot the result (you need matplotlib for this
plt.plot(v,p,'o')
plt.show()

