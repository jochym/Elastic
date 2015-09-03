.. Elastic documentation master file, created by
   sphinx-quickstart on Sun Aug  7 17:19:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Calculation of elastic properties of crystals
=============================================

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.18759.svg
   :target: http://dx.doi.org/10.5281/zenodo.18759

Elastic is a set of python routines for calculation of elastic properties of 
crystals (elastic constants, equation of state, sound velocities, etc.). 
It is a third version of the in-house code I have 
written over several years and is implemented as a extension to the
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ system.
The code was a basis for some of my publications and was 
described briefly in these papers. The code was available to anyone, presented 
at our `Workshop on ab initio Calculations in Geosciences <http://wolf.ifj.edu.pl/workshop/work2008/>`_ 
and used by some of my co-workers but was never properly published with
full documentation, project page etc. Nevertheless the old code is still available
to anyone as `Elastic 2 <http://wolf.ifj.edu.pl/~jochym/elastic2/elastic2.tgz>`_.
I just do not recommend to use it without my help - which I am happy to provide.

In 2010, I have decided to re-implement elastic as a module for the 
`ASE <https://wiki.fysik.dtu.dk/ase/>`_ system and publish it properly under 
the GPL. 

The source code started live on the 
`launchpad project page <https://launchpad.net/elastic>`_ and later in 2014 
moved to the `github repository <https://github.com/jochym/Elastic>`_ with 
corresponding `elastic web page <https://jochym.github.io/Elastic/>`_ and
on-line documentation placed at `Elastic website <http://wolf.ifj.edu.pl/elastic/>`_ 
(you are probably reading from it already). You can obtain the `documentation as a 
PDF file <http://wolf.ifj.edu.pl/~jochym/Elastic.pdf>`_ as well.

The project is open and I welcome patches, ideas and other feedback. 
You can also support the project and motivate me to work on it even more 
by donating using bitcoin address: 1Geq8khANDueVt1QdCS5GU2oNCtdc1uSMv .




.. toctree::
   :maxdepth: 3
   
   intro
   modules
   usage
   endmatter

* :ref:`search`



