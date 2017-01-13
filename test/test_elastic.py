from __future__ import print_function, division
import unittest
import hypothesis
from hypothesis import given, example, assume
from hypothesis.strategies import integers, floats, builds

import elastic

def build_crystal(a,b,c,al,be,ga):
    from ase.atoms import Atoms
    if (al+be+ga > 350) : return
    if (al+be < 1.1*ga) : return
    if (al+ga < 1.1*be) : return
    if (be+ga < 1.1*al) : return
    
    return Atoms('Mg', positions=[(0,0,0)], cell=[a,b,c,al,be,ga],pbc=True)

Crystals = builds(build_crystal,
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

    @given(cr=Crystals)
    def test_lattice_type(self, cr):
        assume(cr is not None)
        assert cr.get_lattice_type() in (1,2,3,4,5,6,7)
        
    @given(cr=Crystals,
            lo=floats(min_value=0.1, max_value=1.0),
            hi=floats(min_value=1.0, max_value=2.0),
            n=integers(min_value=2, max_value=20),
            )
    def test_volume_scan(self, cr, lo, hi, n):
        assume(cr is not None)
        cl=cr.scan_volumes(lo,hi,n,True)
        eps=1e-6
        assert len(cl)==n
        assert abs(cl[0].get_volume() - cr.get_volume()*lo)<eps
        assert abs(cl[-1].get_volume() - cr.get_volume()*hi)<eps
        for c in cl[1:-1]:
            assert cr.get_volume()*lo <= c.get_volume() <= cr.get_volume()*hi


if __name__ == '__main__':
     unittest.main()
