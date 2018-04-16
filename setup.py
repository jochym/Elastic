# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from setuptools_scm import get_version

ver = get_version().split('+')[0].split('.')

if ver[-1].startswith('dev'):
    ver[-1] = ver[-1][3:]
    ver[-2] = str(int(ver[-2])-1)

ver = '.'.join(ver)

print("Building version:", ver)

setup(
    name='elastic',
    version=ver,
    packages=find_packages(),
    license='GPLv3',
    description='Extension for ASE to calculate elastic constants',
    long_description=open('README.rst', 'rb').read().decode('utf-8'),
    long_description_content_type='text/x-rst',
    author='Paweł T. Jochym',
    author_email='Pawel.Jochym@ifj.edu.pl',
    url='https://github.com/jochym/Elastic',
    keywords=['science', 'physics', 'ase', 'elastic constants', 'crystals'],
    requires=['click', 'spglib', 'numpy', 'scipy', 'ase'],
    setup_requires=['docutils', 'sphinx', 'setuptools_scm'],
    provides=['elastic', 'parcalc'],
    platforms=['all'],
    classifiers=[],
    include_package_data=True,
    install_requires=['click', ],
    entry_points='''
        [console_scripts]
        elastic=elastic.cli.elastic:cli
        '''
)
