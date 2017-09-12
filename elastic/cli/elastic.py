#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 1998-2015 by Paweł T. Jochym <pawel.jochym@ifj.edu.pl>
#
#    This file is part of Elastic.

#    Elastic is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    Elastic is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Elastic.  If not, see <http://www.gnu.org/licenses/>.

'''
The elastic command is a command-line tool exposing the functionality
of elastic library for direct use - without writing any python code.
'''

import click
import ase.io
from elastic import get_BM_EOS, get_elastic_tensor
import pkg_resources

verbose = 0


def banner():
    if verbose > 0:
        print('Elastic ver.',
              pkg_resources.get_distribution("elastic").version)


def set_verbosity(v):
    global verbose
    verbose = v


@click.command()
@click.option('-v', '--verbose', count=True)
def log(verbose):
    click.echo('Verbosity: %s' % verbose)
    set_verbosity(verbose)


@click.command()
def cli():
    '''Command-line interface to the elastic library.'''

    banner()


if __name__ == '__main__':
    cli()
