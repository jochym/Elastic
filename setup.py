# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

def scm_config():
    from setuptools_scm.version import (
        postrelease_version,
        get_local_node_and_date,
    )
    return dict(
        version_scheme=postrelease_version,
        local_scheme=get_local_node_and_date,
    )

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
