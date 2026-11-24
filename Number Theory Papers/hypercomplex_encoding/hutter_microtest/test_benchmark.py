import unittest
from benchmark import hamilton,feature_ids,run

class TestHamilton(unittest.TestCase):
    def test_identity(self):
        self.assertEqual(hamilton((1,0,0,0),(7,8,9,10)),(7,8,9,10))
    def test_i_squared(self):
        self.assertEqual(hamilton((0,1,0,0),(0,1,0,0)),(-1,0,0,0))
    def test_ij(self):
        self.assertEqual(hamilton((0,1,0,0),(0,0,1,0)),(0,0,0,1))
    def test_budget(self):
        keys=feature_ids(bytes(range(8)))
        self.assertEqual(len(keys),4)
        self.assertTrue(all(0<=key<512 for key in keys))
    def test_determinism(self):
        data=(b'<page><title>Alpha beta</title><text>Sample text.</text></page>\n')*200
        a=run(data,1000,1200,1200)
        b=run(data,1000,1200,1200)
        self.assertEqual(a['methods'],b['methods'])
        self.assertEqual(a['baseline_ideal_bytes'],b['baseline_ideal_bytes'])

if __name__=='__main__': unittest.main()
