# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("__conda_version__.h") as w:
    for line in w:
        version_t=line.strip()

if None in version:
    print("Failed to get version number in setup.py.")
    raise


setup(
    name='elastic',
    version=version_t,
    packages=find_packages(),
    license='GPLv3',
    description = 'Extension for ASE to calculate elastic constants',
    author = 'Paweł T. Jochym',
    author_email = 'Pawel.Jochym@ifj.edu.pl',
    url = 'https://github.com/jochym/Elastic',
    keywords = ['science', 'physics', 'ase', 'elastic constants', 'crystals'],
    requires = ['spglib','numpy','scipy','ase','docutils','sphinx'],
    provides = ['elastic','parcalc'],
    platforms = ['all'],
    classifiers = [],
)
