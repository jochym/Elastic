Elastic
=======

Elastic is a set of python routines for calculation of elastic properties of 
crystals (elastic constants, equation of state, sound velocities, etc.). 
It is a third version of the in-house code I have 
written over few years and is implemented as a extension to the
[ASE](https://wiki.fysik.dtu.dk/ase/) system.
The code was a basis for some of my publications and was 
described briefly in these papers. The code was available to anyone, presented 
at the [Workshop on ab initio Calculations in Geosciences](http://wolf.ifj.edu.pl/workshop/work2008/)
and used by some of my co-workers. Nevertheless it was never properly published with
full documentation, project page etc. The old code is still available
to anyone as [Elastic 2](http://wolf.ifj.edu.pl/~jochym/elastic2/elastic2.tgz>).
I just do not recommend to use it without my help - which I am happy to provide.

Recently, thanks to generous support of:
* Department of Computational Material Science, Institute of Nuclear Physics, PAN, Poland
* Department of Engineering, University of Saskatchewan, Canada, 

I was able re-implemented elastic as a module for the 
ASE, integrate interface to the QuantumEspresso code and publish 
it properly as a PyPi package (elastic).

The on-line documentation is placed on [Elastic website](http://wolf.ifj.edu.pl/elastic/)
You can obtain the 
[documentation as a PDF file](http://wolf.ifj.edu.pl/~jochym/Elastic.pdf) 
as well.

The installation is simple (provided you have all the dependencies like 
scipy, numpy, ASE, pyspglib to name just a few):

    pip install elastic

If you need a more elaborate 
[installation instruction](http://nbviewer.ipython.org/gist/jochym/a7f552e8b1fced1bc996) 
for computing environment to work with ASE - I have written just such a document.
It is a first in the series of tutorials of this subject and you can 
find it under nbviewer link above.

The project is open and I welcome patches, ideas and other feedback. 

You can also support the project and motivate me to work on it even more 
by donating using bitcoin address: 1Geq8khANDueVt1QdCS5GU2oNCtdc1uSMv .

