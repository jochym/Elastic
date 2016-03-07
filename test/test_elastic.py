from __future__ import print_function
import unittest

class TestElastic(unittest.TestCase):

    def setUp(self):
	    pass
    
    def tearDown(self):
        pass
    
    def test_elastic_import(self):
	    import elastic
	    print('Imported:',elastic.__package__)

    def test_parcalc_import(self):
	    import parcalc
	    print('Imported:',parcalc.__package__)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestElastic)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
