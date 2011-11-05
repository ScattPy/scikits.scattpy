import unittest

from scipy.integrate import quad
from scipy.special import *

from spheroidal import *


class testSpheroidalProlateNorm(unittest.TestCase):
    def test_norm5(self):
        c = 5
        for n in range(0,10):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def test_norm1(self):
        c = 1
        for n in range(0,3):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def test_norm25(self):
        c = 25
        for n in range(0,10):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def test_norm25(self):
        c = 0.2
        for n in range(0,3):
            for m in range(0,n):
                self.check_norm(m, n, c)
        
    def check_norm(self, m, n, c):
        norm = get_pro_norm(m,n,c)
        func = lambda x: power(pro_ang1(m, n, c, x)[0], 2)
        self.assertAlmostEqual((quad(func, -1, 1) / (norm * norm))[0], 1, places = 5)

class testSpheroidalOblateNorm(unittest.TestCase):
    def test_norm5(self):
        c = 5
        for n in range(0,10):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def test_norm1(self):
        c = 1
        for n in range(0,3):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def test_norm25(self):
        c = 25
        for n in range(0,10):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def test_norm25(self):
        c = 0.2
        for n in range(0,3):
            for m in range(0,n):
                self.check_norm(m, n, c)

    def check_norm(self, m, n, c):
        norm = get_obl_norm(m,n,c)
        func = lambda x: power(obl_ang1(m, n, c, x)[0], 2)
        self.assertAlmostEqual((quad(func, -1, 1) / (norm * norm))[0], 1, places = 5)

if __name__ == '__main__':
    unittest.main()
