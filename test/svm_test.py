import unittest

from spheroidal import *
from spheroidal_particles import *
from properties import *

#main integration tests for non-absorbing particle
class testSVM(unittest.TestCase):

    def __init__(self,methodName):
        unittest.TestCase.__init__(self,methodName)
        self.m = 1.3; self.a_b =1.4; self.type = 1; self.xl = 4
        self.particle = Spheroid(self.m,self.a_b,self.type)
        self.particle.set_xl(self.xl)

    #according to (100)
    #holds if k=1 so d should be 2 * c1
    def test_relation1(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_B(3) * svm.get_A(c1,c1,1).T - svm.get_A(c1,c1,3) * svm.get_B(1).T
        right_part = 1j*mat(eye(nmax))
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],right_part[i,k])

    def test_relation2(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_B(3) * svm.get_A(c1,c1,3).T - svm.get_A(c1,c1,3) * svm.get_B(3).T
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],0)

    def test_relation3(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_B(1) * svm.get_A(c1,c1,1).T - svm.get_A(c1,c1,1) * svm.get_B(1).T
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],0)

    def test_relation4(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_A(c1,c1,1).T * svm.get_B(3) - svm.get_A(c1,c1,3).T * svm.get_B(1)
        right_part = 1j*mat(eye(nmax))
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k],right_part[i,k])

    def test_relation5(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_A(c1,c1,3).T * svm.get_A(c1,c1,1) - svm.get_A(c1,c1,1).T * svm.get_A(c1,c1,3)
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k], 0)

    def test_relation6(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_B(3).T * svm.get_B(1) - svm.get_B(1).T * svm.get_B(3)
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k], 0)

    def test_relation7(self):
        particle = self.particle
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        c1 = particle.c1
        relation = svm.get_A(c1,c1,1) * svm.get_B(3).T - svm.get_B(1).T * svm.get_A(c1,c1,3)
        right_part = 1j*mat(eye(nmax))
        for i in range(0,nmax):
            for k in range(0,nmax):
                self.assertAlmostEqual(relation[i,k], right_part[i,k])


    def test_IntegrationCase1(self):
        alpha = pi / 4
        particle = self.particle
        particle.set_xl(1)
        print "c1=" + str(particle.c1)
        print "c2=" + str(particle.c2)
        print "d=" + str(particle.d)
        print "eps="+str(particle.eps)
        print "psi="+str(particle.psi)
        nmax = 8
        svm = SpheroidalSVM(particle,nmax)
        b_sca = svm.getSolution(TMInputWave(alpha))[0]
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
        alpha = pi / 4
        particle = self.particle
        particle.set_xl(3.0)
        print " c1=" + str(particle.c1) + " c2=" + str(particle.c2) + " eps=" + str(particle.eps) + " psi=" + str(particle.psi) + " d="+str(particle.d)
        nmax = 6
        svm = SpheroidalSVM(particle,nmax)
        b_sca = svm.getSolution(TMInputWave(alpha))[0]
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


  