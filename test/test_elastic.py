import unittest

class TestElastic(unittest.TestCase):

    def setUp(self):
	    pass
    
    def tearDown(self):
        pass
    
    def test_basic_import(self):
	    import elastic
	    import parcalc

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestElastic)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
