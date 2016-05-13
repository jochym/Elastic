# -*- coding: utf-8 -*-

#from sys import argv

#if 'upload_docs' in argv :
#    print('upload_docs is not supported by this package.')
#else :
if True :
    from setuptools import setup, find_packages

    with open("__version__.txt") as w:
        for line in w:
            version = line.strip()
            break

    setup(
        name='elastic',
        version=version, 
        packages=find_packages(exclude=['docs']),
        license='GPLv3',
        description = 'Extension for ASE to calculate elastic constants',
        author = 'Paweł T. Jochym',
        author_email = 'Pawel.Jochym@ifj.edu.pl',
        url = 'https://github.com/jochym/Elastic',
        keywords = ['science', 'physics', 'ase', 'elastic constants', 'crystals'],
        requires = ['spglib','numpy','scipy','ase'],
        setup_requires = ['docutils','sphinx'],
        provides = ['elastic','parcalc'],
        platforms = ['all'],
        classifiers = [],
    )
