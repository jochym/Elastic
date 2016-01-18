# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(
    name='elastic',
    use_scm_version=True, 
    packages=find_packages(),
    license='GPLv3',
    description = 'Extension for ASE to calculate elastic constants',
    author = 'Paweł T. Jochym',
    author_email = 'Pawel.Jochym@ifj.edu.pl',
    url = 'https://github.com/jochym/Elastic',
    keywords = ['science', 'physics', 'ase', 'elastic constants', 'crystals'],
    requires = ['spglib','numpy','scipy','ase','docutils','sphinx'],
    setup_requires = ['setuptools_scm'],
    provides = ['elastic','parcalc'],
    platforms = ['all'],
    classifiers = [],
)
