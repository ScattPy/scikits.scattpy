import unittest

from numpy import transpose

from spheroidal import *
from spheroidal_particles import *
from properties import *

#main integration tests for non-absorbing particle
class testSVM(unittest.TestCase):

    #according to (100)
    #holds if k=1 so d should be 2 * c1
    def test_relation1(self):
        nmax = 6
        c1 = 1
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_B(3) * svm.get_A(c1,c1,1).T - svm.get_A(c1,c1,3) * svm.get_B(1).T
        right_part = 1j*mat(eye(nmax))
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],right_part[i,k])

    def test_relation2(self):
        nmax = 6
        c1 = 1.3
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_B(3) * svm.get_A(c1,c1,3).T - svm.get_A(c1,c1,3) * svm.get_B(3).T
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],0)

    def test_relation3(self):
        nmax = 6
        c1 = 1.3
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_B(1) * svm.get_A(c1,c1,1).T - svm.get_A(c1,c1,1) * svm.get_B(1).T
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],0)

    def test_relation4(self):
        nmax = 6
        c1 = 1
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_A(c1,c1,1).T * svm.get_B(3) - svm.get_A(c1,c1,3).T * svm.get_B(1)
        right_part = 1j*mat(eye(nmax))
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],right_part[i,k])

    def test_relation5(self):
        nmax = 6
        c1 = 1.3
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_A(c1,c1,3).T * svm.get_A(c1,c1,1) - svm.get_A(c1,c1,1).T * svm.get_A(c1,c1,3)
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k], 0)

    def test_relation6(self):
        nmax = 6
        c1 = 1.3
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_B(3).T * svm.get_B(1) - svm.get_B(1).T * svm.get_B(3)
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k], 0)

    def test_relation7(self):
        nmax = 6
        c1 = 1
        c2 = 2
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=1)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        relation = svm.get_A(c1,c1,1) * svm.get_B(3).T - svm.get_B(1).T * svm.get_A(c1,c1,3)
        right_part = 1j*mat(eye(nmax))
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k], right_part[i,k])


    def test_IntegrationCase1(self):
        nmax = 6
        c1 = 1
        c2 = 2
        alpha = pi / 3
        k = 1
        particle = ProlateSpheroid(psi=2,c=c1,derivative=0,eps=2)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        b_sca = svm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
        print C_ext, C_sca
        delta = (C_ext-C_sca)/(C_ext+C_sca)
        print delta
        self.assertAlmostEqual(delta,0,5)
        self.assertAlmostEqual(C_ext,C_sca,5)
        self.assertTrue(C_ext>0.)
        self.assertTrue(C_sca>0.)

    def test_IntegrationCase2(self):
        nmax = 10
        c1 = 1
        c2 = 6
        alpha = pi / 3
        k = 1
        particle = OblateSpheroid(psi=2,c=c1,derivative=0,eps=2)
        svm = SpheroidalSVM(particle,c2,c1,nmax)
        b_sca = svm.getSolution(TMInputWave(alpha))[0]
        C_ext = getCext(particle, alpha, k, b_sca, nmax)[0]
        C_sca = getCsca(k, b_sca, nmax)[0]
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


  