# -*- coding: utf-8 -*-

from distutils.core import setup

setup(
    name='elastic',
    packages=['elastic'],
    license='GPLv3',
    description = 'Extension for ASE to calculate elastic constants',
    author = 'Paweł T. Jochym',
    author_email = 'Pawel.Jochym@ifj.edu.pl',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    url = 'https://github.com/jochym/Elastic',
    keywords = ['science', 'physics', 'ase', 'elastic constants', 'crystals'],
    requires = ['pyspglib','numpy','scipy','ase'],
    provides = ['elastic'],
    platforms = ['all'],
    classifiers = [],
)
