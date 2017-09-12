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
import elastic
import pkg_resources
from click import echo

verbose = 0


def banner():
    if verbose > 0:
        echo('Elastic ver. %s\n----------------------' %
             pkg_resources.get_distribution("elastic").version)


def set_verbosity(v):
    global verbose
    verbose = v


def process_calc(fn):
    from time import sleep
    sleep(1)


@click.group()
@click.option('--vasp', 'format', flag_value='vasp',
              help='Use VASP formats (default)', default=True)
@click.option('--abinit', 'format', flag_value='abinit',
              help='Use AbInit formats')
@click.option('--cij', 'action', flag_value='cij',
              help='Generate deformations for Cij (default)', default=True)
@click.option('--eos', 'action', flag_value='eos',
              help='Generate deformations for Equation of State')
@click.option('-v', '--verbose', count=True, help='Increase verbosity')
@click.version_option()
@click.pass_context
def cli(ctx, format, action, verbose):
    '''Command-line interface to the elastic library.'''

    if verbose:
        set_verbosity(verbose)
    banner()


@cli.command()
@click.argument('struct', type=click.Path(exists=True))
@click.pass_context
def gen(ctx, struct):
    '''Generate deformed structures'''

    format = ctx.parent.params['format']
    action = ctx.parent.params['action']
    cryst = ase.io.read(struct, format=format)
    fn_tmpl = action
    if format == 'vasp':
        fn_tmpl += '_%03d.POSCAR'
        kwargs = {'vasp5': True, 'direct': True}
    elif format == 'abinit':
        fn_tmpl += '_%03d.abinit'
        kwargs = {}

    if verbose:
        from elastic.elastic import get_lattice_type
        nr, brav, sg, sgn = get_lattice_type(cryst)
        echo('%s lattice (%s): %s' % (brav, sg, cryst.get_chemical_formula()))

    if action == 'cij':
        systems = elastic.get_elastic_tensor(cryst)
    elif action == 'eos':
        systems = elastic.get_BM_EOS(cryst)

    if verbose:
        echo('Writing %d deformation files.' % len(systems))
    for n, s in enumerate(systems):
        ase.io.write(fn_tmpl % n, s, format=format, **kwargs)


@cli.command()
@click.argument('files', nargs=-1)
@click.pass_context
def proc(ctx, files):
    '''Process calculated structures'''

    format = ctx.parent.params['format']
    action = ctx.parent.params['action']
    with click.progressbar(files) as bar:
        for calc in bar:
            process_calc(calc)


if __name__ == '__main__':
    cli()
