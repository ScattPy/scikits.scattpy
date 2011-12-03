import unittest

from scipy.integrate import quad
from scipy.special import *

from spheroidal import *

class testSpheroidalProlateNorm(unittest.TestCase):
    def test_norm5(self):
        c = 5
        for n in range(1,12):
            for m in range(1,n+1):
                self.check_norm(m, n, c)
            m = 0

    #def test_norm_minus5(self):
    #    c = -5
    #    for n in range(1,4):
    #        for m in range(1,n+1):
    #            self.check_norm(m, n, c)
    #        m = 0
                    
    def test_norm1(self):
        c = 1
        for n in range(1,5):
            for m in range(1,n+1):
                self.check_norm(m, n, c)
            m = 0

    def test_norm20(self):
        self.check_norm(3, 20, 10)

    def test_norm02(self):
        c = 0.2
        for n in range(1,5):
            for m in range(1,n+1):
                self.check_norm(m, n, c)
            m = 0

    def check_norm(self, m, n, c):
        norm = get_pro_norm(m,n,c)
        func = lambda x: power(pro_ang1(m, n, c, x)[0] / norm, 2)
        self.assertAlmostEqual((quad(func, -1, 1, limit = 300))[0], 1, places = 5)

class testSpheroidalOblateNorm(testSpheroidalProlateNorm):
    def check_norm(self, m, n, c):
        norm = get_obl_norm(m,n,c)
        func = lambda x: power(obl_ang1(m, n, c, x)[0] / norm, 2)
        self.assertAlmostEqual((quad(func, -1, 1, limit = 100))[0], 1, places = 5)

class testFactorial(unittest.TestCase):
    
    def test_factorial(self):
        self.assertEqual(get_norm_factorial(1,0), 1)
        self.assertEqual(get_norm_factorial(5,0), 1)
        self.assertEqual(get_norm_factorial(0,0), 1)

        self.assertEqual(get_norm_factorial(3,1), 20)
        self.assertEqual(get_norm_factorial(2,2), 360)


if __name__ == '__main__':
    unittest.main()
