Installation
============

Right now the installation is quite primitive. I will prepare a proper 
packaging for a final release. For now just head to the 
`launchpad code repository <http://bazaar.launchpad.net/~jochym/elastic/trunk/files>`_
and grab two base files: parcalc.py and elastic.py . Then, put them either on 
your python path or in the working directory. There are no special installation
procedures required. You just have to have fully functional ASE system with a
working quantum mechanical calculator. I am trying to write the code 
independently from any particular setup (i.e. my own setup).

To use the code you need to install following python packages (most, if not all, are available in major linux distributions):

* `SciPy and NumPy <http://www.scipy.org/>`_ libraries
* `matplotlib <http://matplotlib.sourceforge.net/>`_ (not strictly required,
  but needed for testing and plotting)
* `ASE <https://wiki.fysik.dtu.dk/ase/>`_ system
* Some ASE calculator (VASP, GPAW, abinit, ...), but be warned that for now 
  the code was developed using VASP only. I will be happy to help you extending
  it to other calculators.
* `spglib <http://spglib.sourceforge.net/>`_ space group library 
* `pyspglib <http://spglib.sourceforge.net/pyspglibForASE/>`_ python space group module

Testing
-------

All modules have small testing sets at the end. You can run these test by 
simply running each module as a python script::

    python parcalc.py

which will run a short series of single-point calculations on the MgO unit
cell and plot the resulting equation of state. The main module testing routine::

    python elastic.py

will run the equation of state and elastic tensor calculations for a set of 
small crystals - one for each Bravais lattice. This may take some considerable
time. The testing routines will probably not work out of the box in your system.
Review the comments at the end of the files to make them work. I will try to make 
them as setup-agnostic as possible.

Usage
=====

Once you have everything installed and running you can run your first real 
calculation. The testing code at the end of the elastic.py may be used as 
an example how to do it. The first step is to import the modules to your 
program (the examples here use VASP calculator)::

    from ase.lattice.spacegroup import crystal
    from parcalc import ClusterVasp
    import ase.units as units

next we need to create the crystal, MgO in this case::

    a = 4.194
    MgO = crystal(['Mg', 'O'], 
                    [(0, 0, 0), (0.5, 0.5, 0.5)], 
                    spacegroup=225,
                    cellpar=[a, a, a, 90, 90, 90])

We need a calculator for our job, here we use VASP and ClusterVasp defined 
in the parcalc module. You can probably replace this calculator by any other ASE
calculator but this was not tested yet. Thus let us define the calculator::

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
                isif=4, 
                nsw=30,
                ediff=1e-8, 
                ibrion=2)
    calc.set(kpts=[3,3,3])
    
    # Set the calculation mode first.
    # Full structure optimization in this case.
    # Not all calculators have this type of internal minimizer!
    calc.set(isif=3)

Finally, run our first calculation. Obtain relaxed structure and 
residual pressure after optimization::

    print "Residual pressure: %.3f GPa" % (
                cryst.get_pressure(cryst.get_stress())/units.GPa)


