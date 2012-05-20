import unittest

from spheroidal_particles import *
from properties import *

#main integration tests for non-absorbing particle
from spheroidal_pmm import SpheroidalPMM

class testPMM(unittest.TestCase):

    def test_IntegrationCase1(self):
        alpha = pi / 4
        m=1.5;a_b=1.2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1)
        nmax = 6
        pmm = SpheroidalPMM(particle,nmax)
        b_sca = pmm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        print C_ext, C_sca
        delta = (C_ext-C_sca)/(C_ext+C_sca)
        print delta
        self.assertAlmostEqual(delta,0,5)
        self.assertAlmostEqual(C_ext,C_sca,5)
        self.assertTrue(C_ext>0.)
        self.assertTrue(C_sca>0.)

    def test_IntegrationCase2(self):
        alpha = pi / 3
        m=1.3;a_b=1.1;type=-1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1)
        nmax = 6
        pmm = SpheroidalPMM(particle,nmax)
        b_sca = pmm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        print C_ext, C_sca
        delta = (C_ext-C_sca)/(C_ext+C_sca)
        print delta
        self.assertAlmostEqual(delta,0,5)
        self.assertAlmostEqual(C_ext,C_sca,5)
        self.assertTrue(C_ext>0.)
        self.assertTrue(C_sca>0.)

    def test_IntegrationCase3(self):
        alpha = pi / 3
        m=1.3;a_b=1.2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1)
        nmax = 6
        pmm = SpheroidalPMM(particle,nmax)
        b_sca = pmm.getMatrixSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        print C_ext, C_sca
        delta = (C_ext-C_sca)/(C_ext+C_sca)
        print delta
        self.assertAlmostEqual(delta,0,5)
        self.assertAlmostEqual(C_ext,C_sca,5)
        self.assertTrue(C_ext>0.)
        self.assertTrue(C_sca>0.)

    def test_IntegrationCase4(self):
        alpha = pi / 3
        m=1.3;a_b=1.2;type=-1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1)
        nmax = 6
        pmm = SpheroidalPMM(particle,nmax)
        b_sca = pmm.getMatrixSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, b_sca, nmax)[0]
        C_sca = getCsca(particle, b_sca, nmax)[0]
        print C_ext, C_sca
        delta = (C_ext-C_sca)/(C_ext+C_sca)
        print delta
        self.assertAlmostEqual(delta,0,5)
        self.assertAlmostEqual(C_ext,C_sca,5)
        self.assertTrue(C_ext>0.)
        self.assertTrue(C_sca>0.)

if __name__ == '__main__':
    #import cProfile
    #cProfile.run('unittest.main()','profiler_results')
    unittest.main()


  