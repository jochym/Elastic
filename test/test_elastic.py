from __future__ import print_function
import unittest

class TestElastic(unittest.TestCase):

    def setUp(self):
	    pass
    
    def tearDown(self):
        pass
    
    def test_elastic_import(self):
	    import elastic

    def test_parcalc_import(self):
	    import parcalc

if __name__ == '__main__':
     unittest.main()
