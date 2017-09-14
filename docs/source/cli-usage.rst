Installation
============

.. image:: https://anaconda.org/jochym/elastic/badges/version.svg
   :target: https://anaconda.org/jochym/elastic
.. image:: https://anaconda.org/jochym/elastic/badges/installer/conda.svg
   :target: https://conda.anaconda.org/jochym

Conda
-------

.. highlight:: bash

The installation procedure is quite simple if you use, *highly recommended*
`conda package manager <http://conda.pydata.org/miniconda.html>`_::

    conda install -c conda-forge elastic

The above command installs elastic with all dependencies into your current
conda environment. If you want to add my anaconda.org channel into your conda
installation you need to run following command::

    conda config --add channels conda-forge

The above method has additional benefit of providing current installation of
ASE and spglib libraries.

PyPi
------

The package is published simultaneously on conda and pypi. The second 
recommended way to install elastic is with pip::

    pip install elastic

which should install the package and all its dependencies. 
Note that the number of dependencies is rather large and some of them are
fairly large. Most of them, however, are just standard scientific python
packages - almost all are present in standard anaconda install.


Manual
--------

To install the code *pedestrian way* you need to install following python 
packages (most, if not all, are available in major linux distributions):

* `SciPy and NumPy <http://www.scipy.org/>`_ libraries
* `matplotlib <http://matplotlib.org/>`_ (not strictly required,
  but needed for testing and plotting)
* `ASE <https://wiki.fysik.dtu.dk/ase/>`_ system
* `spglib <https://atztogo.github.io/spglib/>`_ space group library 
* Some ASE calculator (VASP, GPAW, abinit, ...), but be warned 
  that for now the code was developed using VASP only. I will be happy to 
  help you extending it to other calculators. The code can be used without
  supported ASE calculator using command line interface and external, 
  independent calculation tool.

This is highly system-dependent and I am unable to provide detailed support for
this type of install - I use conda install of ASE/elastic myself!

Some legacy `installation guides <https://github.com/jochym/qe-doc/blob/master/Installation.ipynb>`_ 
which may help you with manual process could be find at the 
`QE-doc project pages <https://jochym.github.io/qe-doc/>`_.

Essentially, you need to clone the repository and run::

    python setup.py install

in the main directory. But this is really not recommended way to install.
Use it at your own risk and if you know what you are doing.

Testing
-------

The simplest verification if everything works is just using the ``elastic``
utility to see the help screen::

    elastic --help

or version of the package::

    elastic --version

Additionally the whole package has a set of unittests based on hypothesis
package. The tests are self-contained and do not require any external packages
like DFT programs (e.g. VASP). You can run these tests by executing
following command in the source directory::

    python -m unittest discover -s test -b


Usage
=====

Starting from ver. 5.0, the command line utility ``elastic`` is a primary
interface to the package and the direct python programming with the library
is relegated to non-standard calculations. 
The basic calculation scheme can be summarized with the following list:

* Prepare the basic structure data in the form of VASP POSCAR file or 
  abinit input file. Support for other programs can be added relatively
  easily. Contact the author if you need it. 
  The structure should be fully optimized represent what you consider to
  be ground state of the system.

* run ``elastic`` on the structure to generate deformed structures probing
  the properties of the system::

    elastic -v --cij gen -n 5 -s 2 POSCAR

  which generates a set of deformed systems named ``cij_XXX.POSCAR``, where
  ``XXX`` is replaced by numbers with 000 corresponding to undisturbed structure.

* run your DFT program (VASP, abinit, etc.) on all systems. This step depends
  on your particular system, and you need to handle it yourself. You need to 
  make sure that for each system the internal degrees of freedom are 
  optimized and full stress tensor gets calculated. Example 
  bash script handling this task on my cluster looks like this::

    #!/bin/bash
    
    # Command to run vasp in current directory
    RUN_VASP="/opt/bin/run-vasp541"

    for s in $* ; do
        d=${s%%.POSCAR}
        echo -n  $d ": "
        mkdir calc-$d
        (
            cd calc-$d
            ln -s ../INCAR ../KPOINTS ../POTCAR ../vasprun.conf .
            ln -s ../$s POSCAR
            $RUN_VASP
        )
    done

  This produces a set of directories: ``calc-cij_XXX`` with completed 
  single-point calculations.

* run ``elastic`` again to post-process the calculations. 
  We do that by feeding it with output from the DFT calculations.
  Remember to put undisturbed structure at the first position::

    elastic -v --cij proc calc-cij_000/vasprun.xml calc-cij_*/vasprun.xml

.. highlight:: console

You can test this procedure using data provided as a reference in the 
``tests/data`` directory. If you run the script on the provided data you
should get following output::

    elastic -v --cij proc calc-cij_000/vasprun.xml calc-cij_*/vasprun.xml
    
    Cij solution
    ------------------------------
     Solution rank:  3
     Square of residuals: 0.00053
     Relative singular values:
     1.0000   0.7071   0.6354  

    Elastic tensor (GPa):
       C_11     C_12     C_44  
    ------------------------------
     321.15    95.88   143.44 

The data provided correspond to cubic MgO crystal. The DFT calculation setup
is tuned to provide quick results for testing and *should not* be used as
a guide for production calculations. You *need* to determine proper
calculation setup for your system.
