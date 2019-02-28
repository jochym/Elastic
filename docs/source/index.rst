.. Elastic documentation master file, created by
   sphinx-quickstart on Sun Aug  7 17:19:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Calculation of elastic properties of crystals
=============================================

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.18759.svg
   :target: http://dx.doi.org/10.5281/zenodo.18759

------

**New version 5.0 released**

The new version is API *incompatible* with the previous versions.
It provides a new command line utility as the main user interface
to the package - which hoppefully will broaden the user base byond
python users.

------

Elastic is a set of python routines for calculation of elastic properties of 
crystals (elastic constants, equation of state, sound velocities, etc.). 
It is a fifth version of the in-house code I have 
written over several years and is implemented as a extension to the
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ system and a script providing
interface to the library not requiring knowledge of python or ASE system.
The code was a basis for some of my publications and was 
described briefly in these papers. The code was available to anyone, presented 
at our `Workshop on ab initio Calculations in Geosciences <http://wolf.ifj.edu.pl/workshop/work2008/>`_ 
and used by some of my co-workers but was never properly published with
full documentation, project page etc. In 2010, I have decided to re-implement 
elastic as a module for the `ASE <https://wiki.fysik.dtu.dk/ase/>`_ system 
and publish it properly under the GPL as versions 3 and 4. 

Later, in 2017, needs of users nudged me into implementing the command-line 
front-end to the library which is included with version 5.0 of the package.

The version 5.0 also changes API of the library from mix-in class to the set
of simple functions providing functionality of the module. The workflow of the 
package was also changed to prepare data - calculate - post-process style, 
which is better suited to serious research work.

The source code started live on the 
`launchpad project page <https://launchpad.net/elastic>`_ and later in 2014 
moved to the `github repository <https://github.com/jochym/Elastic>`_ with 
corresponding `elastic web page <https://jochym.github.io/Elastic/>`_ and
on-line documentation placed at `Elastic website <http://wolf.ifj.edu.pl/elastic/>`_ 
(you are probably reading from it already).

The project is free software and I welcome patches, ideas and other feedback.


.. toctree::
   :maxdepth: 3
   
   intro
   cli-usage
   lib-usage
   modules
   endmatter

* :ref:`search`



