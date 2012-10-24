
import unittest
import cowpy

class TestReduce(unittest.TestCase):
    def setUp(self):
        self.domain = cowpy.DistributedDomain([16,16])
        print "running test on %d processes" % self.domain.cart_size

    def tearDown(self):
        del self.domain

    def test_sum(self):
        val = 1
        res = self.domain.reduce(val, type=int, op='sum')
        self.assertEqual(type(res), int)
        self.assertEqual(res, self.domain.cart_size)

    def test_max(self):
        val = self.domain.cart_rank
        res = self.domain.reduce(val, type=float, op='max')
        self.assertEqual(type(res), float)
        self.assertEqual(res, self.domain.cart_size - 1)


if __name__ == "__main__":
    unittest.main(exit=False)
