#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright 2011 by Pawel T. Jochym <pawel.jochym@ifj.edu.pl>
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

.. _par-calc-mod:

Parallel Calculator Module
^^^^^^^^^^^^^^^^^^^^^^^^^^

Parallel calculator module is an extension of the standard 
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ calculator working in the
parallel cluster environment. It is very useful in all situations where 
you need to run several, independent calculations and you have a large 
cluster of machines at your disposal (probably with some queuing system).

This implementation uses VASP but the code can be easily adapted for use
with other ASE calculators with minor changes.
The final goal is to provide a universal module for parallel 
calculator execution in the cluster environment.

The SIESTA code by Georgios Tritsaris <gtritsaris@seas.harvard.edu>
Not fully tested after merge.

'''

from .parcalc import ClusterVasp, ClusterSiesta, ParCalculate


