Elastic |DOI|
=============

|Build Status| |Binstar Badge| |Binstar Badge| |Binstar Badge| |Binstar
Badge| |Research software impact|

Elastic is a set of python routines for calculation of elastic
properties of crystals (elastic constants, equation of state, sound
velocities, etc.). It is a third version of the in-house code I have
written over few years and is implemented as a extension to the
`ASE <https://wiki.fysik.dtu.dk/ase/>`__ system. The code was a basis
for some of my publications and was described briefly in these papers.
The code was available to anyone, presented at the `Workshop on ab
initio Calculations in
Geosciences <http://wolf.ifj.edu.pl/workshop/work2008/>`__ and used by
some of my co-workers. Nevertheless it was never properly published with
full documentation, project page etc. The old code is still available to
anyone as `Elastic
2 <http://wolf.ifj.edu.pl/~jochym/elastic2/elastic2.tgz%3E>`__. I just
do not recommend to use it without my help - which I am happy to
provide.

Recently, thanks to generous support of: \* Department of Computational
Material Science, Institute of Nuclear Physics, PAN, Poland \*
Department of Engineering, University of Saskatchewan, Canada,

I was able re-implemented elastic as a module for the ASE, integrate
interface to the QuantumEspresso code and publish it properly as a PyPi
package (elastic).

The on-line documentation is placed on
`ReadTheDocs <http://elastic.rtfd.org/>`__ or `Elastic
website <http://wolf.ifj.edu.pl/elastic/>`__. You can obtain the
`documentation as a PDF
file <https://media.readthedocs.org/pdf/elastic/stable/elastic.pdf>`__
as well.

Installation |Binstar Badge|
----------------------------

The installation is simple if you are using conda package menager:

::

    conda install -c jochym elastic

If you need a more elaborate `installation
instruction <http://nbviewer.ipython.org/github/jochym/qe-doc/blob/master/Installation.ipynb>`__
for computing environment to work with ASE - I have written just such a
document. It is a first in the series of tutorials of this subject and
you can find it under nbviewer link above.

The project is open and I welcome patches, ideas and other feedback.

You can also support the project and motivate me to work on it even more
by donating using bitcoin address: 1Geq8khANDueVt1QdCS5GU2oNCtdc1uSMv .

.. |DOI| image:: https://zenodo.org/badge/doi/10.5281/zenodo.18759.svg
   :target: http://dx.doi.org/10.5281/zenodo.18759
.. |Build Status| image:: https://travis-ci.org/jochym/Elastic.svg?branch=master
   :target: https://travis-ci.org/jochym/Elastic
.. |Binstar Badge| image:: https://anaconda.org/jochym/elastic/badges/version.svg
   :target: https://anaconda.org/jochym/elastic
.. |Binstar Badge| image:: https://anaconda.org/jochym/elastic/badges/build.svg
   :target: https://anaconda.org/jochym/elastic/builds
.. |Binstar Badge| image:: https://anaconda.org/jochym/elastic/badges/downloads.svg
   :target: https://anaconda.org/jochym/elastic
.. |Binstar Badge| image:: https://anaconda.org/jochym/elastic/badges/license.svg
   :target: https://anaconda.org/jochym/elastic
.. |Research software impact| image:: http://depsy.org/api/package/pypi/elastic/badge.svg
   :target: http://depsy.org/package/python/elastic
.. |Binstar Badge| image:: https://anaconda.org/jochym/elastic/badges/installer/conda.svg
   :target: https://conda.anaconda.org/jochym
