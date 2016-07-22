# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from setuptools_scm import get_version

ver = get_version().split('.')
if len(ver)>3 :
    # Post release version:
    #  - add patchlevel at the end
    #  - decrease back last number
    pl=int(ver[3].split('+')[0][3:])
    ver[2]=str(int(ver[2])-1)
    ver=ver[:3]
    ver.append(str(pl))

ver='.'.join(ver)

setup(
    name='elastic',
    version=ver,
    packages=find_packages(),
    license='GPLv3',
    description = 'Extension for ASE to calculate elastic constants',
    author = 'Paweł T. Jochym',
    author_email = 'Pawel.Jochym@ifj.edu.pl',
    url = 'https://github.com/jochym/Elastic',
    keywords = ['science', 'physics', 'ase', 'elastic constants', 'crystals'],
    requires = ['click','spglib','numpy','scipy','ase'],
    setup_requires = ['docutils','sphinx','setuptools_scm'],
    provides = ['elastic','parcalc'],
    platforms = ['all'],
    classifiers = [],
)
