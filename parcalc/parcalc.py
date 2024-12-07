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

from __future__ import print_function, division
import logging
import ase
from ase.calculators.vasp import Vasp
from ase.calculators.siesta import Siesta
from ase.calculators.aims import Aims
from ase.calculators.calculator import Calculator, FileIOCalculator, all_changes

try :  # Python3
    from queue import Empty
except ImportError : # Python2
    from Queue import Empty

from multiprocessing import Process, Queue

import time
import os
import tempfile
import shutil
from copy import deepcopy
from subprocess import check_output
from contextlib import contextmanager

class _NonBlockingRunException(Exception):
    '''
    Internal exception. Should never be propagated outside.
    '''
    def __str__(self):
        return '''The __NonBlockingRunException should be caught inside
        the calculator class. If you got it outside it is a bug.
        Contact the author and/or submit a bug ticket at github.'''

from traceback import print_stack

@contextmanager
def work_dir(path):
    '''
    Context menager for executing commands in some working directory.
    Returns to the previous wd when finished.

    Usage:
    >>> with work_dir(path):
    ...   subprocess.call('git status')

    '''
    starting_directory = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(starting_directory)


class ClusterVasp(Vasp):
    '''
    Adaptation of VASP calculator to the cluster environment where you often
    have to make some preparations before job submission. You can easily
    adapt this class to your particular environment. It is also easy to
    use this as a template for other type of calculator.
    '''

    def __init__(self, nodes=1, ppn=8, block=True, ncl=False, **kwargs):
        Vasp.__init__(self, **kwargs)
        self.nodes=nodes
        self.ppn=ppn
        self.block=block
        self.ncl=ncl
        self.calc_running=False
        self.working_dir=os.getcwd()

    def prepare_calc_dir(self):
        '''
        Prepare the calculation directory for VASP execution.
        This needs to be re-implemented for each local setup.
        The following code reflects just my particular setup.
        '''
        with open("vasprun.conf","w") as f:
            f.write('NODES="nodes=%s:ppn=%d"\n' % (self.nodes, self.ppn))
            f.write('BLOCK=%d\n' % (self.block,))
            if self.ncl :
                f.write('NCL=%d\n' % (1,))
            #print(self.nodes, self.ppn)

    def calc_finished(self):
        '''
        Check if the lockfile is in the calculation directory.
        It is removed by the script at the end regardless of the
        success of the calculation. This is totally tied to
        implementation and you need to implement your own scheme!
        '''
        #print_stack(limit=5)
        if not self.calc_running :
            #print('Calc running:',self.calc_running)
            return True
        else:
            # The calc is marked as running check if this is still true
            # We do it by external scripts. You need to write these
            # scripts for your own system.
            # See examples/scripts directory for examples.
            with work_dir(self.working_dir) :
                o=check_output(['check-job'])
            #print('Status',o)
            if o[0] in b'R' :
                # Still running - we do nothing to preserve the state
                return False
            else :
                # The job is not running maybe it finished maybe crashed
                # We hope for the best at this point ad pass to the
                # Standard update function
                return True

    def set(self,**kwargs):
        if 'block' in kwargs :
            self.block=kwargs['block']
            del kwargs['block']
        else :
            self.block=True

        if 'ncl' in kwargs :
            self.ncl=kwargs['ncl']
            del kwargs['ncl']
        else :
            self.ncl=False
        Vasp.set(self, **kwargs)

    def clean(self):
        with work_dir(self.working_dir) :
            Vasp.clean(self)

    def update(self, atoms):
        if self.calc_running :
            # we have started the calculation and have
            # nothing to read really. But we need to check
            # first if this is still true.
            if  self.calc_finished():
                # We were running but recently finished => read the results
                # This is a piece of copy-and-paste programming
                # This is a copy of code from Vasp.calculate
                self.calc_running=False
                with work_dir(self.working_dir) :
                    atoms_sorted = ase.io.read('CONTCAR', format='vasp')
                    if self.int_params['ibrion'] > -1 and self.int_params['nsw'] > 0:
                        # Update atomic positions and unit cell with the ones read
                        # from CONTCAR.
                        atoms.positions = atoms_sorted[self.resort].positions
                        atoms.cell = atoms_sorted.cell
                    self.converged = self.read_convergence()
                    Vasp.set_results(self,atoms)
                return
            else :
                return
        # We are not in the middle of calculation.
        # Update as normal
        Vasp.update(self, atoms)

    def set_results(self, atoms):
        with work_dir(self.working_dir) :
            #print('set_results')
            Vasp.set_results(self, atoms)

    def run(self):
        '''
        Blocking/Non-blocing run method.
        In blocking mode it just runs parent run method.
        In non-blocking mode it raises the __NonBlockingRunException
        to bail out of the processing of standard calculate method
        (or any other method in fact) and signal that the data is not
        ready to be collected.
        '''
        # This is only called from self.calculate - thus
        # we do not need to change to working_dir
        # since calculate already did
        Vasp.run(self)
        if not self.block :
            #print('Interrupt processing of calculate', os.getcwd())
            raise _NonBlockingRunException

    def calculate(self, atoms):
        '''
        Blocking/Non-blocking calculate method

        If we are in blocking mode we just run, wait for
        the job to end and read in the results. Easy ...

        The non-blocking mode is a little tricky.
        We need to start the job and guard against it reading
        back possible old data from the directory - the queuing
        system may not even started the job when we get control
        back from the starting script. Thus anything we read
        after invocation is potentially garbage - even if it
        is a converged calculation data.

        We handle it by custom run function above which
        raises an exception after submitting the job.
        This skips post-run processing in the calculator, preserves
        the state of the data and signals here that we need to wait
        for results.
        '''

        with work_dir(self.working_dir) :
            self.prepare_calc_dir()
            self.calc_running=True
            #print('Run VASP.calculate')
            try :
                Vasp.calculate(self, atoms)
                self.calc_running=False
                #print('VASP.calculate returned')
            except _NonBlockingRunException as e:
                # We have nothing else to docs
                # until the job finishes
                #print('Interrupted ', self.working_dir, os.getcwd())
                pass



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
        with open("siestarun.conf","w") as f:
            f.write('NODES="nodes=%d:ppn=%d"' % (self.nodes, self.ppn))
            #print(self.nodes, self.ppn)

    def get_potential_energy(self, atoms):
        self.prepare_calc_dir()
        Siesta.get_potential_energy(self, atoms)

    def clean(self):
        self.converged = False
        return

class ClusterAims(Aims):
    '''
    Encapsulating Aims calculator for the cluster environment.
    '''
    def __init__(self, nodes=1, ppn=8, **kwargs):
        Aims.__init__(self, **kwargs)
        self.nodes=nodes
        self.ppn=ppn

    def prepare_calc_dir(self):
        with open("siestarun.conf","w") as f:
            f.write('NODES="nodes=%d:ppn=%d"' % (self.nodes, self.ppn))
            #print(self.nodes, self.ppn)

    def run(self):
        self.prepare_calc_dir()
        Aims.run(self)


class RemoteCalculator(Calculator):
    '''
    Remote calculator based on ASE calculator class.
    This class is only involved with the machanics of remotly executing
    the software and transporting the data. The calculation is
    delegated to the actual calculator class.
    '''
    
    # Queue system submit command
    qsub_tool='qsub'
    qstat_tool='qstat'
    
    qsub_cmd='cd %(rdir)s ; %(qsub_tool)s -N %(title)s -l procs=%(procs)d ./run-pw.pbs'

    # Remote execution command
    remote_exec_cmd='ssh %(user)s@%(host)s "%(command)s"'

    # If you cannot mount the data directory into your system it is best
    # to use the rsync command to transfer the results back into the system.

    # Command for copying the data out to the computing system
    copy_out_cmd='rsync -a "%(ldir)s" "%(user)s@%(host)s:%(rdir)s"'
    # Command for copying the data in after the calculation
    copy_in_cmd='rsync -a "%(user)s@%(host)s:%(rdir)s" "%(ldir)s"'

    # Template for the PBS batch job
    pbs_template=''

    # Command to check the state of the job
    pbs_check_cmd='''%(qstat_tool)s -f %(jobid)s |grep job_state |awk '{print $3}' '''
    
    # Access data
    host=''
    user=''
    
    # Location:
    # local working directory
    wdir='.'
    # Remote working directory relative to the home directory or absolute
    rdir='.'

    # Repetition timer (seconds) for checkin the state of the job.
    job_check_time=15
    
    def __init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                 atoms=None, calc=None, block=False, **kwargs):
        '''Basic calculator implementation.
        restart: str
            Prefix for restart file.  May contain a directory.  Default
            is None: don't restart.
        ignore_bad_restart_file: bool
            Ignore broken or missing restart file.  By default, it is an
            error if the restart file is missing or broken.
        label: str
            Name used for all files.  May contain a directory.
        atoms: Atoms object
            Optional Atoms object to which the calculator will be
            attached.  When restarting, atoms will get its positions and
            unit-cell updated from file.

        Create a remote execution calculator based on actual ASE calculator
        calc.
        '''
        logging.debug("Calc: %s Label: %s" % (calc, label))
        Calculator.__init__(self, restart, ignore_bad_restart_file, label, atoms, **kwargs)
        logging.debug("Dir: %s Ext: %s" % (self.directory, self.ext))
        self.calc=calc
        self.jobid=None
        self.block=block

    def write_pbs_in(self,properties):
        with work_dir(self.directory):
            with open(os.path.join(self.directory,'run-ase-calc.pbs'),'w') as fh:
                fh.write(self.pbs_template % {
                    'command': self.build_command(self,prop=properties,
                                                    params=self.parameters)
                    })

    def build_command(self, prop=None, params=None):
        '''
        Build command strings for local or remote execution
        '''
        if prop is None :
            prop = ['energy']
        if params is None :
            params = {}
        cmd=self.qsub_cmd % {
                'qsub_tool': self.qsub_tool,
                'qstat_tool': self.qstat_tool,
                'title': self.label,
                'procs': self.parameters['procs'],
                'rdir': os.path.join(self.parameters['rdir'],
                                     os.path.split(self.directory)[-1])
            }
        cmd=self.remote_exec_cmd % {
                'command': cmd,
                'user': self.parameters['user'],
                'host': self.parameters['host']
            }
        return cmd

    def write_input(self, atoms=None, properties=None, system_changes=all_changes):
        '''Write input file(s).'''
        if properties is None :
            properties = ['energy']
        with work_dir(self.directory):
            self.calc.write_input(self, atoms, properties, system_changes)
            self.write_pbs_in(properties)
            subprocess.call(self.copy_out_cmd % {
                                'ldir': self.directory,
                                'rdir': self.parameters['rdir'],
                                'user': self.parameters['user'],
                                'host': self.parameters['host']
                            }, shell=True)

    def job_ready(self):
        try :
            cmd=self.remote_exec_cmd % {
                'command': self.pbs_check_cmd % {
                                    'qsub_tool': self.qsub_tool,
                                    'qstat_tool': self.qstat_tool,
                                    'jobid':self.jobid
                                    },
                'user': self.parameters['user'],
                'host': self.parameters['host']
                }
            state=subprocess.check_output(cmd, shell=True).split()[-1]
        except (subprocess.CalledProcessError, IndexError) :
            # Unknown state. We assume it has finished and continue
            state='N'

        return not (state in ['Q','R'])


    def run_calculation(self, atoms=None, properties=None,
                            system_changes=all_changes):
        '''
        Internal calculation executor. We cannot use FileIOCalculator
        directly since we need to support remote execution.
        
        This calculator is different from others.
        It prepares the directory, launches the remote process and
        raises the exception to signal that we need to come back for results
        when the job is finished.
        '''
        if properties is None :
            properties = ['energy']
        self.calc.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        if self.command is None:
            raise RuntimeError('Please configure Remote calculator!')
        olddir = os.getcwd()
        errorcode=0
        try:
            os.chdir(self.directory)
            output = subprocess.check_output(self.command, shell=True)
            self.jobid=output.split()[0]
            self.submited=True
            #print "Job %s submitted. Waiting for it." % (self.jobid)
            # Waiting loop. To be removed.
        except subprocess.CalledProcessError as e:
            errorcode=e.returncode
        finally:
            os.chdir(olddir)
        
        if errorcode:
            raise RuntimeError('%s returned an error: %d' %
                               (self.name, errorcode))
        self.read_results()


    def read_results(self):
        """Read energy, forces, ... from output file(s)."""
        
        if self.submited:
            # The job has been submitted. Check the state.
            if not self.job_ready() :
                if self.block :
                    while not self.job_ready() :
                        time.sleep(self.job_check_time)
                else :
                    raise CalcNotReadyError
                    

            # Assume the calc finished. Copy the files back.
            subprocess.call(self.copy_in_cmd % {
                'ldir': self.wdir,
                'rdir': os.path.join(self.parameters['rdir'],os.path.split(self.directory)[-1]),
                'user': self.parameters['user'],
                'host': self.parameters['host']
            }, shell=True)
                
        
        fn=os.path.join(self.directory,'pw.out')
        # Read the pan-ultimate line of the output file
        try:
            ln=open(fn).readlines()[-2]
            if ln.find('JOB DONE.')>-1 :
                # Job is done we can read the output
                r=read_quantumespresso_textoutput(fn)
                self.submited=False
                self.jobid=None
            else :
                # Job not ready.
                raise CalcNotReadyError
        except (IOError, IndexError) :
            # Job not ready.
            raise CalcNotReadyError

        # All is fine - really read the results
        self.calc.read_results(self)

    @classmethod
    def ParallelCalculate(cls,syslst,properties=None,system_changes=all_changes):
        '''
        Run a series of calculations in parallel using (implicitely) some
        remote machine/cluster. The function returns the list of systems ready
        for the extraction of calculated properties.
        '''
        print('Launching:',end=' ')
        sys.stdout.flush()
        if properties is None :
            properties = ['energy']
        for n,s in enumerate(syslst):
            try :
                s.calc.block=False
                s.calc.calculate(atoms=s,properties=properties,system_changes=system_changes)
            except CalcNotReadyError:
                s.calc.block=True
            print(n+1, end=' ')
            sys.stdout.flush()
        print()
        print('     Done:', end=' ')
        sys.stdout.flush()
        for n,s in enumerate(syslst):
            s.calc.read_results()
            print( n+1, end=' ')
            sys.stdout.flush()
        print()
        return syslst



verbose=True


class __PCalcProc(Process):
    '''
    Internal helper class representing the calculation process isolated
    from the rest of the ASE script. The process (not thread) runs in
    the separate directory, created on-the-fly and removed at the end
    if the cleanup is true and we are in blocking (default) mode.
    In this mode it is vital for the calculator to read in all the
    results after the run since the files will be removed as soon as the
    "calculate" function terminates. You can pass False to the cleanup
    argument to prevent the clean-up. This is very usefull for debuging.

    The process waits without any time-out for the calculation to finish.
    This is great for short and simple calculations and quick testing.
    You run the job and get back your results.
    '''

    def __init__(self, iq, oq, calc, prefix, cleanup=True):
        Process.__init__(self)
        self.calc=calc
        self.basedir=os.getcwd()
        self.place=tempfile.mkdtemp(prefix=prefix, dir=self.basedir)
        self.iq=iq
        self.oq=oq
        self.CleanUp=cleanup

    def run(self):
        with work_dir(self.place) :
            n,system=self.iq.get()
            system.set_calculator(deepcopy(self.calc))
            system.get_calculator().block=True
            system.get_calculator().working_dir=self.place
            #print("Start at :", self.place)
            if hasattr(self.calc, 'name') and self.calc.name=='Siesta':
                system.get_potential_energy()
            else:
                system.get_calculator().calculate(system)

            #print("Finito: ", os.getcwd(), system.get_volume(), system.get_pressure())
            self.oq.put([n,system])
            if self.CleanUp :
                system.get_calculator().clean()
                os.chdir(self.basedir)
                shutil.rmtree(self.place, ignore_errors=True)


def ParCalculate(systems,calc,cleanup=True,block=True,prefix="Calc_"):
    '''
    Run calculators in parallel for all systems.
    Calculators are executed in isolated processes and directories.
    The resulting objects are returned in the list (one per input system).
    '''

    if isinstance(systems, list) :
        sysl = systems
    else :
        sysl = [systems]

    if block :
        iq=Queue(len(sysl)+1)
        oq=Queue(len(sysl)+1)

        # Create workers
        for s in sysl:
            __PCalcProc(iq, oq, calc, prefix=prefix, cleanup=cleanup).start()

        # Put jobs into the queue
        for n,s in enumerate(sysl):
            iq.put([n,s])
            # Protection against too quick insertion
            time.sleep(0.2)

        if verbose :
            print("Workers started:", len(sysl))

       # Collect the results
        res=[]
        while len(res)<len(sysl) :
            n,s=oq.get()
            res.append([n,s])
            #print("Got from oq:", n, s.get_volume(), s.get_pressure())
    else :
        # We do not need the multiprocessing complications for non-blocking
        # workers. We just run all in sequence.
        basedir=os.getcwd()
        res=[]
        for n,s in enumerate(sysl):
            s.set_calculator(deepcopy(calc))
            s.get_calculator().block=block
            place=tempfile.mkdtemp(prefix=prefix, dir=basedir)
            os.chdir(place)
            s.get_calculator().working_dir=place
            #print("Start at :", place)
            if hasattr(calc, 'name') and calc.name=='Siesta':
                s.get_potential_energy()
            else:
                s.get_calculator().calculate(s)
            os.chdir(basedir)
            #print("Submited", s.get_calculator().calc_finished(), os.getcwd())
            # Protection against too quick insertion
            time.sleep(0.2)
            res.append([n,s])
        if verbose :
            print("Workers started:", len(sysl))

    return [r for ns,s in enumerate(sysl) for nr,r in res if nr==ns]



# Testing routines using VASP as a calculator in the cluster environment.
# TODO: Make it calculator/environment agnostic
if __name__ == '__main__':
    from ase.lattice.spacegroup import crystal
    from ase.units import GPa
    import elastic
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

    print("Residual pressure: %.3f GPa" % (MgO.get_pressure()/GPa))
    calc.clean()

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

