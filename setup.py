# -*- coding: utf-8 -*-

from distutils.core import setup

setup(
    name='elastic',
    version='4.0.5',
    packages=['elastic'],
    license='GPLv3',
    description = 'Extension for ASE to calculate elastic constants',
    author = 'Paweł T. Jochym',
    author_email = 'Pawel.Jochym@ifj.edu.pl',
    url = 'https://github.com/jochym/Elastic',
    download_url = 'https://github.com/jochym/Elastic/tarball/4.0.2', 
    keywords = ['science', 'physics', 'ase', 'elastic constants', 'crystals'],
    requires = ['pyspglib','numpy','scipy','ase'],
    provides = ['elastic'],
    platforms = ['all'],
    classifiers = [],
)
