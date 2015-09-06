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

'''

from __future__ import print_function, division, absolute_import


from ase.atoms import Atoms
from .elastic import BMEOS, Crystal, ElasticCrystal


# Enhance the Atoms class by adding new capabilities
for k in ElasticCrystal.__dict__ :
    if k[:2]!='__' and k[-2:]!='__' :
        #print('Implanting', k)
        setattr(Atoms, k, ElasticCrystal.__dict__[k])
#    Atoms.__atoms_init__=Atoms.__init__
#    Atoms.__init__=__Crystal.__crystal_init__



