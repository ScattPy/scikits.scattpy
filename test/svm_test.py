import unittest

from spheroidal import *
from spheroidal_particles import *
from properties import *

#main integration tests for non-absorbing particle
class testSVM(unittest.TestCase):
    def test_svm1(self):
        nmax = 5
        c1 = 1.0
        c2 = 1.0
        alpha = pi / 4
        particle = ProlateSpheroid()
        b_sca = getSolution(particle, TMInputWave(alpha), c1, c2, nmax)[0]
        k = 1
        C_ext = getCext(particle, alpha, 1, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        self.assertAlmostEqual(C_ext, C_sca, 7)

    def test_svm2(self):
        nmax = 3
        c1 = 2.0
        c2 = 2.0
        alpha = pi / 2
        particle = OblateSpheroid()
        b_sca = getSolution(particle, TMInputWave(alpha), c1, c2, nmax)[0]
        k = 4
        C_ext = getCext(particle, alpha, 1, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        self.assertAlmostEqual(C_ext, C_sca, 7)

if __name__ == '__main__':
    unittest.main()


  