Implementation
==============

Elastic is implemented as an extension module to ASE system

The Elastic package provides, basically, one main python module and one 
auxiliary module ( :ref:`par-calc-mod` ) which can be useful outside of 
the scope of the main code. The :ref:`par-calc-mod` is not distributed 
separately but can be just placed by itself somewhere in your python path
and used with any part of the `ASE <https://wiki.fysik.dtu.dk/ase/>`_ .
I hope it will be incorporated in the main project sometime in the future.

Additionally, the package provides an executable program for those that 
do not wish to write their own python code. This program provides basic
standard procedure for elastic tensor calculation. It is described in the
separate section and in the examples (not yet included in the repository).


.. _modules:

Modules
-------

.. _par-calc-mod:

Parallel Calculator Module
^^^^^^^^^^^^^^^^^^^^^^^^^^

Parallel calculator module is an extension of the standard 
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ calculator working in the
parallel cluster environment. It is very useful in all situations where 
you need to run several, independent calculations and you have a large 
cluster of machines at your disposal (probably with some queuing system).

.. automodule:: parcalc
   :members:


.. _elastic-mod:

Elastic Module
^^^^^^^^^^^^^^

This module depends on :ref:`par-calc-mod` for parallelisation of 
independent calculations.

.. automodule:: elastic
   :members:



