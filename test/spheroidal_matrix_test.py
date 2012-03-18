import unittest
from spheroidal import *

#tests without using derivatives
#prolate spheroid
from spheroidal_particles import TMInputWave, Spheroid
from spheroidal_svm import *

class testSpheroidalMatrixA(unittest.TestCase):

    def test_getA11Simplest(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 3
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_A11()
        cv = get_cv(1,1,c,particle.type)
        self.assertAlmostEqual(A[0,0],rad3_cv(1,1,c,cv,particle.type,psi)[0])
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1,1],rad3_cv(1,2,c,cv,particle.type,psi)[0])
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2,2],rad3_cv(1,3,c,cv,particle.type,psi)[0])
        self.assertAlmostEqual(A[1,0], 0)
        self.assertAlmostEqual(A[0,1], 0)
        self.assertAlmostEqual(A[0,2], 0)
        self.assertAlmostEqual(A[2,0], 0)
        self.assertAlmostEqual(A[2,1], 0)
        self.assertAlmostEqual(A[1,2], 0)

    def test_getA12Simplest(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 3
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_A12()
        cv = get_cv(1,1,c,particle.type)
        self.assertAlmostEqual(A[0,0],-rad1_cv(1,1,c,cv,particle.type,psi)[0])
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1,1],-rad1_cv(1,2,c,cv,particle.type,psi)[0])
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2,2],-rad1_cv(1,3,c,cv,particle.type,psi)[0])
        self.assertAlmostEqual(A[1,0], 0)
        self.assertAlmostEqual(A[0,1], 0)
        self.assertAlmostEqual(A[0,2], 0)
        self.assertAlmostEqual(A[2,0], 0)
        self.assertAlmostEqual(A[2,1], 0)
        self.assertAlmostEqual(A[1,2], 0)

    def test_getA21Simplest(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 3
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_A21()
        cv = get_cv(1,1,c,particle.type)
        coeff = lambda nu: metric_phi(nu,particle) * metric_nu(nu,particle) / metric_psi(nu,particle) * ang1_cv(1, 1, c, cv, particle.type, nu)[0] * ang1_cv(1, 1, c, cv, particle.type, nu)[0]
        coeff = quad(coeff,-1,1)
        self.assertAlmostEqual(A[0,0],rad3_cv(1,1,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1,1],rad3_cv(1,2,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2,2],rad3_cv(1,3,c,cv,particle.type,psi)[1] * coeff)
        self.assertAlmostEqual(A[1,0], 0)
        self.assertAlmostEqual(A[0,1], 0)
        self.assertAlmostEqual(A[0,2], 0)
        self.assertAlmostEqual(A[2,0], 0)
        self.assertAlmostEqual(A[2,1], 0)
        self.assertAlmostEqual(A[1,2], 0)

    def test_getA22Simplest(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 3
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_A22()
        cv = get_cv(1,1,c,particle.type)
        coeff = lambda nu: metric_phi(nu,particle) * metric_nu(nu,particle) / metric_psi(nu,particle) * ang1_cv(1, 1, c, cv, particle.type, nu)[0] * ang1_cv(1, 1, c, cv, particle.type, nu)[0]
        coeff = - quad(coeff,-1,1)
        self.assertAlmostEqual(A[0,0],rad1_cv(1,1,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1,1],rad1_cv(1,2,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2,2],rad1_cv(1,3,c,cv,particle.type,psi)[1] * coeff)
        self.assertAlmostEqual(A[1,0], 0)
        self.assertAlmostEqual(A[0,1], 0)
        self.assertAlmostEqual(A[0,2], 0)
        self.assertAlmostEqual(A[2,0], 0)
        self.assertAlmostEqual(A[2,1], 0)
        self.assertAlmostEqual(A[1,2], 0)

    def test_getA10Simplest(self):
        m=2;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 3
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_A10()
        cv = get_cv(1,1,c,particle.type)
        self.assertAlmostEqual(A[0,0],-rad1_cv(1,1,c,cv,particle.type,psi)[0])
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1,1],-rad1_cv(1,2,c,cv,particle.type,psi)[0])
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2,2],-rad1_cv(1,3,c,cv,particle.type,psi)[0])
        self.assertAlmostEqual(A[1,0], 0)
        self.assertAlmostEqual(A[0,1], 0)
        self.assertAlmostEqual(A[0,2], 0)
        self.assertAlmostEqual(A[2,0], 0)
        self.assertAlmostEqual(A[2,1], 0)
        self.assertAlmostEqual(A[1,2], 0)

    def test_getA20Simplest(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 3
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_A20()
        cv = get_cv(1,1,c,particle.type)
        coeff = lambda nu: metric_phi(nu,particle) * metric_nu(nu,particle) / metric_psi(nu,particle) * ang1_cv(1, 1, c, cv, particle.type, nu)[0] * ang1_cv(1, 1, c, cv, particle.type, nu)[0]
        coeff = - quad(coeff,-1,1)
        self.assertAlmostEqual(A[0,0],rad1_cv(1,1,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1,1],rad1_cv(1,2,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2,2],rad1_cv(1,3,c,cv,particle.type,psi)[1] * coeff)
        self.assertAlmostEqual(A[1,0], 0)
        self.assertAlmostEqual(A[0,1], 0)
        self.assertAlmostEqual(A[0,2], 0)
        self.assertAlmostEqual(A[2,0], 0)
        self.assertAlmostEqual(A[2,1], 0)
        self.assertAlmostEqual(A[1,2], 0)

    def test_get_fullA1(self):
        m=2;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 1
        svm = SpheroidalSVM(particle,nmax)
        A = svm.get_fullA()
        self.assertAlmostEqual(A[0, 0], svm.get_A11())
        self.assertAlmostEqual(A[0, 1], svm.get_A12())
        self.assertAlmostEqual(A[1, 0], svm.get_A21())
        self.assertAlmostEqual(A[1, 1], svm.get_A22())

    def test_get_fullA2(self):
        m=2;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 2
        svm = SpheroidalSVM(particle,nmax)
        c = particle.c1
        psi = particle.psi
        A = svm.get_fullA()
        A11 = svm.get_A11()
        A12 = svm.get_A12()
        A21 = svm.get_A21()
        A22 = svm.get_A22()
        self.assertAlmostEqual(A[0, 0], A11[0, 0])
        self.assertAlmostEqual(A[0, 1], A11[0, 1])
        self.assertAlmostEqual(A[1, 0], A11[1, 0])
        self.assertAlmostEqual(A[1, 1], A11[1, 1])
        self.assertAlmostEqual(A[0, 2], A12[0, 0])
        self.assertAlmostEqual(A[0, 3], A12[0, 1])
        self.assertAlmostEqual(A[1, 2], A12[1, 0])
        self.assertAlmostEqual(A[1, 3], A12[1, 1])
        self.assertAlmostEqual(A[2, 0], A21[0, 0])
        self.assertAlmostEqual(A[2, 1], A21[0, 1])
        self.assertAlmostEqual(A[3, 0], A21[1, 0])
        self.assertAlmostEqual(A[3, 1], A21[1, 1])
        self.assertAlmostEqual(A[2, 2], A22[0, 0])
        self.assertAlmostEqual(A[2, 3], A22[0, 1])
        self.assertAlmostEqual(A[3, 2], A22[1, 0])
        self.assertAlmostEqual(A[3, 3], A22[1, 1])


class testSpheroidalMatrixZ(unittest.TestCase):

    def test_getFullB1(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 1
        svm = SpheroidalSVM(particle,nmax)
        B = svm.get_fullB()
        self.assertAlmostEqual(B[0], svm.get_A10())
        self.assertAlmostEqual(B[1], svm.get_A20())

    def test_getFullB2(self):
        m=1;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nmax = 2
        svm = SpheroidalSVM(particle,nmax)
        B = svm.get_fullB()
        B1 = svm.get_A10()
        B2 = svm.get_A20()
        self.assertAlmostEqual(B[0, 0], B1[0, 0])
        self.assertAlmostEqual(B[0, 1], B1[0, 1])
        self.assertAlmostEqual(B[1, 0], B1[1, 0])
        self.assertAlmostEqual(B[1, 1], B1[1, 1])
        self.assertAlmostEqual(B[2, 0], B2[0, 0])
        self.assertAlmostEqual(B[2, 1], B2[0, 1])
        self.assertAlmostEqual(B[3, 0], B2[1, 0])
        self.assertAlmostEqual(B[3, 1], B2[1, 1])