Elastic
=======

![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/jochym/Elastic/pypi.yml)
![PyPI - Version](https://img.shields.io/pypi/v/elastic)
![PyPI - Downloads](https://img.shields.io/pypi/dm/elastic?label=PyPi)
![Conda Version](https://img.shields.io/conda/vn/conda-forge/elastic)
![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/elastic?label=forge)
![GitHub License](https://img.shields.io/github/license/jochym/Elastic)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.593721.svg)](https://doi.org/10.5281/zenodo.593721)
[![Documentation Status](https://readthedocs.org/projects/elastic/badge/?version=latest)](https://elastic.readthedocs.io/en/latest/?badge=latest)

Elastic is a set of python routines for calculation of elastic
properties of crystals (elastic constants, equation of state, sound
velocities, etc.). It is a third version of the in-house code I have
written over few years and is implemented as a extension to the
[ASE](https://wiki.fysik.dtu.dk/ase/) system. The code was a basis for
some of my publications and was described briefly in these papers. The
code was available to anyone, presented at the [Workshop on ab initio
Calculations in Geosciences](http://wolf.ifj.edu.pl/workshop/work2008/)
and used by some of my co-workers. This code is a re-implementation of
elastic as a module for the ASE.

The on-line documentation is placed on
[ReadTheDocs](http://elastic.rtfd.org/) or [Elastic
website](http://wolf.ifj.edu.pl/elastic/). You can obtain the
[documentation as a PDF
file](https://media.readthedocs.org/pdf/elastic/stable/elastic.pdf) as
well.


Installation
------------

The installation is simple if you are using conda package menager:

    conda install -c conda-forge elastic

If you need a more elaborate [installation
instruction](http://nbviewer.ipython.org/github/jochym/qe-doc/blob/master/Installation.ipynb)
for computing environment to work with ASE - I have written just such a
document. It is a first in the series of tutorials of this subject and
you can find it under nbviewer link above.

The project is open and I welcome patches, ideas and other feedback.

Acknowledgments
---------------

In the past the project was partially supported by:

-   Department of Computational Material Science, Institute of Nuclear
    Physics, PAN, Poland
-   Department of Engineering, University of Saskatchewan, Canada
