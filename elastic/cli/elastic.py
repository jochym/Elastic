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

from __future__ import print_function, absolute_import, division
import click
import ase.io
import elastic
import pkg_resources
from click import echo

verbose = 0


def banner():
    if verbose > 1:
        echo('Elastic ver. %s\n----------------------' %
             pkg_resources.get_distribution("elastic").version)


def set_verbosity(v):
    global verbose
    verbose = v


def process_calc(fn):
    from time import sleep
    sleep(1)


@click.group()
@click.option('--vasp', 'frmt', flag_value='vasp',
              help='Use VASP formats (default)', default=True)
@click.option('--abinit', 'frmt', flag_value='abinit',
              help='Use AbInit formats')
@click.option('--cij', 'action', flag_value='cij',
              help='Generate deformations for Cij (default)', default=True)
@click.option('--eos', 'action', flag_value='eos',
              help='Generate deformations for Equation of State')
@click.option('-v', '--verbose', count=True, help='Increase verbosity')
@click.version_option()
@click.pass_context
def cli(ctx, frmt, action, verbose):
    '''Command-line interface to the elastic library.'''

    if verbose:
        set_verbosity(verbose)
    banner()


@cli.command()
@click.option('-n', '--num', 'num', default=5, type=int,
              help='Number of generated deformations per axis (default: 5)')
@click.option('-l', '--lo', 'lo', default=0.98, type=float,
              help='Lower relative volume for EOS scan (default: 0.98)')
@click.option('-h', '--hi', 'hi', default=1.02, type=float,
              help='Upper relative volume for EOS scan (default: 1.02)')
@click.option('-s', '--size', 'size', default=2.0, type=float,
              help='Deformation size for Cij scan (% or deg., default: 2.0)')
@click.argument('struct', type=click.Path(exists=True))
@click.pass_context
def gen(ctx, num, lo, hi, size, struct):
    '''Generate deformed structures'''

    frmt = ctx.parent.params['frmt']
    action = ctx.parent.params['action']
    cryst = ase.io.read(struct, format=frmt)
    fn_tmpl = action
    if frmt == 'vasp':
        fn_tmpl += '_%03d.POSCAR'
        kwargs = {'vasp5': True, 'direct': True}
    elif frmt == 'abinit':
        fn_tmpl += '_%03d.abinit'
        kwargs = {}

    if verbose:
        from elastic.elastic import get_lattice_type
        nr, brav, sg, sgn = get_lattice_type(cryst)
        echo('%s lattice (%s): %s' % (brav, sg, cryst.get_chemical_formula()))
        if action == 'cij':
            echo('Generating {:d} deformations of {:.1f}(%/degs.) per axis'.format(
                    num, size))
        elif action == 'eos':
            echo('Generating {:d} deformations from {:.3f} to {:.3f} of V0'.format(
                    num, lo, hi))

    if action == 'cij':
        systems = elastic.get_elementary_deformations(cryst, n=num, d=size)
    elif action == 'eos':
        systems = elastic.scan_volumes(cryst, n=num, lo=lo, hi=hi)

    systems.insert(0, cryst)
    if verbose:
        echo('Writing %d deformation files.' % len(systems))
    for n, s in enumerate(systems):
        ase.io.write(fn_tmpl % n, s, format=frmt, **kwargs)


@cli.command()
@click.argument('files', type=click.Path(exists=True), nargs=-1)
@click.pass_context
def proc(ctx, files):
    '''Process calculated structures'''

    def calc_reader(fn, verb):
        if verb:
            echo('Reading: {:<60s}\r'.format(fn), nl=False, err=True)
        return ase.io.read(fn)

    action = ctx.parent.params['action']
    systems = [calc_reader(calc, verbose) for calc in files]
    if verbose :
        echo('', err=True)
    if action == 'cij':
        cij = elastic.get_elastic_tensor(systems[0], systems=systems[1:])
        msv = cij[1][3].max()
        eps = 1e-4
        if verbose:
            echo('Cij solution\n'+30*'-')
            echo(' Solution rank: {:2d}{}'.format(
                    cij[1][2],
                    ' (undetermined)' if cij[1][2] < len(cij[0]) else ''))
            if cij[1][2] == len(cij[0]):
                echo(' Square of residuals: {:7.2g}'.format(cij[1][1]))
            echo(' Relative singular values:')
            for sv in cij[1][3]/msv:
                echo('{:7.4f}{}'.format(
                        sv, '* ' if (sv) < eps else '  '), nl=False)
            echo('\n\nElastic tensor (GPa):')
            for dsc in elastic.elastic.get_cij_order(systems[0]):
                echo('{: >7s}  '.format(dsc), nl=False)
            echo('\n'+30*'-')
        for c, sv in zip(cij[0], cij[1][3]/msv):
            echo('{:7.2f}{}'.format(
                    c/ase.units.GPa, '* ' if sv < eps else '  '), nl=False)
        echo()
    elif action == 'eos':
        eos = elastic.get_BM_EOS(systems[0], systems=systems[1:])
        eos[1] /= ase.units.GPa
        if verbose:
            echo('# %7s (A^3)  %7s (GPa)   %7s' % ("V0", "B0", "B0'"))
        echo('    %7.2f        %7.2f        %7.2f' % tuple(eos))


if __name__ == '__main__':
    cli()
