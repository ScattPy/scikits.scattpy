import unittest

from spheroidal import *
from spheroidal_particles import *
from properties import *


class testSVM(unittest.TestCase):
    def test_svm(self):
        nmax = 5
        c1 = 1.0
        c2 = 1.0
        alpha = pi/4
        particle = ProlateSpheroid()
        b_sca = getSolution(particle,TMInputWave(pi/4),c1, c2, nmax)[0]

        k = 1
        C_ext = getCext(particle, alpha, 1, b_sca, nmax)[0]
        C_sca = getCsca(k,b_sca,nmax)[0]
        self.assertAlmostEqual(C_ext, C_sca, 3)

if __name__ == '__main__':
    unittest.main()


  