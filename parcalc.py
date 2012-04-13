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

Class description
"""""""""""""""""
'''

from ase.calculators.vasp import *
from ase.calculators.siesta import *
from Queue import Empty
from multiprocessing import Process, Queue

import time
import os
import tempfile
import shutil

class ClusterVasp(Vasp):
    '''
    Adaptation of VASP calculator to the cluster environment where you often
    have to make some preparations before job submission. You can easily 
    adapt this class to your particular environment. It is also easy to
    use this as a template for other type of calculator.
    '''
    
    def __init__(self, nodes=1, ppn=8, **kwargs):
        Vasp.__init__(self, **kwargs)
        self.nodes=nodes
        self.ppn=ppn
        
    def prepare_calc_dir(self):
        '''
        Prepare the calculation directory for VASP execution.
        This needs to be re-implemented for each local setup.
        The following code reflects just my particular setup.
        '''
        f=open("vasprun.conf","w")
        f.write('NODES="nodes=%d:ppn=%d"' % (self.nodes, self.ppn))
        #print  self.nodes, self.ppn
        f.close()
   
    def calculate(self, atoms):
        self.prepare_calc_dir()
        Vasp.calculate(self, atoms)



class ClusterSiesta(Siesta):
    '''
    Siesta calculator. Not fully tested by me - so this should be considered
    beta quality. Nevertheless it is based on working implementation
    '''
    def __init__(self, nodes=1, ppn=8, **kwargs):
        Siesta.__init__(self, **kwargs)
        self.nodes=nodes
        self.ppn=ppn
    
    def prepare_calc_dir(self):
        f=open("siestarun.conf","w")
        f.write('NODES="nodes=%d:ppn=%d"' % (self.nodes, self.ppn))
        #print  self.nodes, self.ppn
        f.close()
    
    def get_potential_energy(self, atoms):
        self.prepare_calc_dir()
        Siesta.get_potential_energy(self, atoms)

    def clean(self):
        Siesta.converged = False
        return

verbose=True

def ParCalculate(systems,calc,cleanup=True,prefix="Calc_"):
    '''
    Run calculators in parallel for all systems. 
    Calculators are executed in isolated processes and directories.
    The resulting objects are returned in the list (one per input system).
    '''

    class PCalcProc(Process):
        '''
        Internal Class representing the calculation process isolated
        from the rest of the ASE script. The process (not thread) runs in 
        the separate, temporary directory, created on-the-fly and removed at
        the end. It is vital for the calculator to read in all the results 
        after the run since the files will be removed as soon as the 
        "calculate" function terminates. You can pass False to the cleanup
        argument to prevent the clean-up. This is very usefull for debuging.
        '''
        
        def __init__(self, iq, oq, calc, prefix, cleanup):
            Process.__init__(self)
            self.calc=calc
            self.basedir=os.getcwd()
            self.place=tempfile.mkdtemp(prefix=prefix, dir=self.basedir)
            self.iq=iq
            self.oq=oq
            self.CleanUp=cleanup
        
        def run(self):
            wd=os.getcwd()
            os.chdir(self.place)
            system=self.iq.get()
            system.set_calculator(self.calc)
            if hasattr(self.calc, 'name') and self.calc.name=='Siesta':
                system.get_potential_energy()
            else:
                self.calc.calculate(system)
            #print "Finito: ", self.place, os.getcwd()
            self.oq.put(system)
            #print system.get_volume(), system.get_isotropic_pressure(system.get_stress())
            if self.CleanUp :
                self.calc.clean()
                os.chdir(wd)
                shutil.rmtree(self.place, ignore_errors=True)


    if type(systems) != type([]) :
        sys=[systems]
    else :
        sys=systems

    iq=Queue(len(sys)+1)
    oq=Queue(len(sys)+1)
        
    # Create workers    
    for s in sys:
        PCalcProc(iq, oq, calc, prefix=prefix, cleanup=cleanup).start()

    # Put jobs into the queue
    for s in sys:
        iq.put(s)
        # Protection against too quick insertion
        time.sleep(0.2)
    
    time.sleep(2)
    if verbose : 
        print len(sys), "Workers started"
    
   # Collect the results
    res=[]
    while len(res)<len(sys) :
        s=oq.get()
        res.append(s)
        #print "Got from oq:", s.get_volume(), s.get_isotropic_pressure(s.get_stress())
    return res

# Testing routines using VASP as a calculator in the cluster environment.
# TODO: Make it calculator/environment agnostic
if __name__ == '__main__':
    from ase.lattice.spacegroup import crystal
    import numpy
    from pylab import *

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
    
    print MgO.get_isotropic_pressure(MgO.get_stress())
                   
    sys=[]
    for av in numpy.linspace(a*0.95,a*1.05,5):
        sys.append(crystal(['Mg', 'O'], [(0, 0, 0), (0.5, 0.5, 0.5)], spacegroup=225,
                   cellpar=[av, av, av, 90, 90, 90]))
                       
    pcalc=ClusterVasp(nodes=1,ppn=8)
    pcalc.set(prec = 'Accurate', xc = 'PBE', lreal = False, isif=2, nsw=20, ibrion=2, kpts=[1,1,1])
    res=ParCalculate(sys,pcalc)
    
    v=[]
    p=[]
    for s in res :
        v.append(s.get_volume())
        p.append(s.get_isotropic_pressure(s.get_stress()))
    
    plot(v,p,'o')
    show()
    
