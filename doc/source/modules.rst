Implementation
==============

Elastic is implemented as an extension module to 
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ system

The Elastic package provides, basically, one main python module and one 
auxiliary module (:ref:`par-calc-mod`) which can be useful outside of 
the scope of the main code. The :ref:`par-calc-mod` is not distributed 
separately but can be just placed by itself somewhere in your python path
and used with any part of the ASE.
I hope it will be incorporated in the main project sometime in the future.

Additionally, the package provides an executable program for those that 
do not wish to write their own python code. This program provides basic
standard procedure for elastic tensor calculation. It is described in the
separate section and in the examples (not yet included in the repository).


.. _modules:

Modules
-------

.. automodule:: parcalc
   :members:



.. automodule:: elastic
   :members:



