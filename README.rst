Elastic
=======

|Build Status| |PVersion Badge| |CVersion Badge| |Downloads Badge| |License Badge| |Research software impact| |DOI| |Codacy Badge|

Elastic is a set of python routines for calculation of elastic
properties of crystals (elastic constants, equation of state, sound
velocities, etc.). It is a third version of the in-house code I have
written over few years and is implemented as a extension to the
`ASE <https://wiki.fysik.dtu.dk/ase/>`__ system. The code was a basis
for some of my publications and was described briefly in these papers.
The code was available to anyone, presented at the 
`Workshop on ab initio Calculations in Geosciences <http://wolf.ifj.edu.pl/workshop/work2008/>`__ 
and used by some of my co-workers. This code is a re-implementation
of elastic as a module for the ASE.

The on-line documentation is placed on
`ReadTheDocs <http://elastic.rtfd.org/>`__ or 
`Elastic website <http://wolf.ifj.edu.pl/elastic/>`__. You can obtain the
`documentation as a PDF file <https://media.readthedocs.org/pdf/elastic/stable/elastic.pdf>`__
as well.

The docs are also published at: 
`Elastic over IPFS <https://ipfs.io/ipns/QmSHUr59SLvWEZq7URTDGgouVaq7vFJYG7HqtBgL4s6M2u/>`__

The stable site hash is: `/ipns/QmSHUr59SLvWEZq7URTDGgouVaq7vFJYG7HqtBgL4s6M2u/`

Installation 
-------------

The installation is simple if you are using conda package menager:

::

    conda install -c conda-forge elastic

If you need a more elaborate 
`installation instruction <http://nbviewer.ipython.org/github/jochym/qe-doc/blob/master/Installation.ipynb>`__
for computing environment to work with ASE - I have written just such a
document. It is a first in the series of tutorials of this subject and
you can find it under nbviewer link above.

The project is open and I welcome patches, ideas and other feedback.

Acknowledgments
---------------

In the past the project was partially supported by:

- Department of Computational Material Science, Institute of Nuclear Physics, PAN, Poland
- Department of Engineering, University of Saskatchewan, Canada

.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.18759.svg
   :target: http://dx.doi.org/10.5281/zenodo.18759
.. |Build Status| image:: https://travis-ci.org/jochym/Elastic.svg?branch=master
   :target: https://travis-ci.org/jochym/Elastic
.. |PVersion Badge| image:: https://img.shields.io/pypi/v/elastic.svg   
   :target: https://pypi.org/project/elastic/
.. |CVersion Badge| image:: https://anaconda.org/conda-forge/elastic/badges/version.svg
   :target: https://anaconda.org/conda-forge/elastic
.. |Downloads Badge| image:: https://anaconda.org/jochym/elastic/badges/downloads.svg
   :target: https://anaconda.org/jochym/elastic
.. |License Badge| image:: https://anaconda.org/jochym/elastic/badges/license.svg
   :target: https://anaconda.org/jochym/elastic
.. |Research software impact| image:: http://depsy.org/api/package/pypi/elastic/badge.svg
   :target: http://depsy.org/package/python/elastic
.. |Codacy Badge| image:: https://api.codacy.com/project/badge/Grade/9f348d724d564ce399beb009ab9a3654
   :alt: Codacy Badge
   :target: https://www.codacy.com/app/jochym/Elastic?utm_source=github.com&utm_medium=referral&utm_content=jochym/Elastic&utm_campaign=badger
