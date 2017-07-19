from __future__ import print_function, division
import unittest
import hypothesis
from hypothesis import given, example, assume
from hypothesis.strategies import integers, floats, builds
from numpy import pi, allclose, cos, sin, array, zeros

import elastic

def ctg(x):
    return cos(x)/sin(x)

def csc(x):
    return 1/sin(x)


def build_crystal(a,b,c,al,be,ga):
    from ase.atoms import Atoms
    if (al+be+ga > 350) : return
    if (al+be < 1.1*ga) : return
    if (al+ga < 1.1*be) : return
    if (be+ga < 1.1*al) : return
    
    atms=Atoms('Mg', positions=[(0,0,0)], cell=[a,b,c,al,be,ga],pbc=True)
    atms._test_data=[a,b,c,al*pi/180,be*pi/180,ga*pi/180]
    return atms

crystals = builds(build_crystal,
                   a=floats(min_value=1.0, max_value=20.0),
                   b=floats(min_value=1.0, max_value=20.0),
                   c=floats(min_value=1.0, max_value=20.0),
                   al=floats(min_value=15.0, max_value=120.0),
                   be=floats(min_value=15.0, max_value=120.0),
                   ga=floats(min_value=15.0, max_value=120.0),
                )


class TestElastic(unittest.TestCase):

    def setUp(self):
        pass
    
    def tearDown(self):
        pass
    
    def test_elastic_import(self):
        import elastic

    def test_parcalc_import(self):
        import parcalc

    def test_mixin_action(self):
        import elastic
        import ase
        for s in elastic.ElasticCrystal.__dict__:
            if not s.startswith('__'):
                assert (s in ase.atoms.Atoms.__dict__)

    @given(cr=crystals)
    def test_lattice_type(self, cr):
        assume(cr is not None)
        assert cr.get_lattice_type() in (1,2,3,4,5,6,7)

    @given(cr=crystals)
    def test_lattice_vecang(self, cr):
        assume(cr is not None)
        assert allclose(cr.get_vecang_cell(),cr._test_data)

    @given(cr=crystals,
            ax=integers(min_value=0,max_value=2),
            prc=floats(min_value=-90,max_value=90),
            )
    def test_deform_shrink(self,cr,ax,prc):
        assume(cr is not None)
        dc=cr.get_deformed_cell(axis=ax, size=prc)
        assert allclose((dc.get_volume()/cr.get_volume()-1)*100 , prc)
        assert allclose((dc.get_cell_lengths_and_angles()[ax]/cr._test_data[ax] - 1)*100, prc)

    @given(cr=crystals,
            ax=integers(min_value=0,max_value=2),
            size=floats(min_value=-10,max_value=10),
            )
    def test_deform_angels(self,cr,ax,size):
        assume(cr is not None)
        delta = zeros(3)
        delta[ax]= pi*size/180
        (alp, bet, gam) = array(cr._test_data[3:]) + delta
        # Guard against non-existing lattice vectors
        t =1 - (ctg(bet)*ctg(gam)-cos(alp)*csc(bet)*csc(gam))**2
        assume(t > 1e-6)
        dc=cr.get_deformed_cell(axis=ax+3, size=size)
        assert allclose(pi*dc.get_cell_lengths_and_angles()[ax+3]/180,
                        cr._test_data[ax+3] + pi*size/180, atol=pi*1e-3/180)

    @given(cr=crystals,
            lo=floats(min_value=0.1, max_value=1.0),
            hi=floats(min_value=1.0, max_value=2.0),
            n=integers(min_value=2, max_value=20),
            )
    def test_volume_scan(self, cr, lo, hi, n):
        assume(cr is not None)
        cl=cr.scan_volumes(lo,hi,n,True)
        eps=1e-6
        assert len(cl)==n
        assert allclose(cl[0].get_volume() , cr.get_volume()*lo)
        assert allclose(cl[-1].get_volume() , cr.get_volume()*hi)
        for c in cl[1:-1]:
            assert cr.get_volume()*lo <= c.get_volume() <= cr.get_volume()*hi


if __name__ == '__main__':
     unittest.main()
