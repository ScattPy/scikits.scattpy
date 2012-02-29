import unittest
from spheroidal import *

#tests without using derivatives
#prolate spheroid
from spheroidal_particles import TMInputWave

class Particle:
    type = 1

    def __init__(self, psi, d=1.5, eps=2):
        self.psi = psi
        self.d = d
        self.eps = eps

    def function(self, nu):
        return self.psi

    def derivative(self, nu):
        return 0


class testSpheroidalMatrixA(unittest.TestCase):
    def test_get_a_functionR1_1(self):
        n = 1;
        m = 0;
        c = 2.5;
        z = 0.5;
        type = 1;
        particle = Particle(2.5)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 1, particle)(z), -0.1280857374177537)

    def test_get_a_functionR1_2(self):
        n = 3;
        m = 2;
        c = 2.5;
        z = -0.5;
        type = 1;
        particle = Particle(1.5)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 1, particle)(z), -0.1547625029193905)

    def test_get_a_functionR1_3(self):
        n = 5;
        m = 0;
        c = 0.5;
        z = -0.9;
        type = 1;
        particle = Particle(5)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 1, particle)(z), 0.000642762007655166)

    def test_get_a_functionR1_4(self):
        n = 4;
        m = 2;
        c = 4.5;
        z = -0.6;
        type = 1;
        psi = 6
        particle = Particle(psi)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 1, particle)(z), 0.03417296826676850, 5)

    def test_get_a_functionR3_1(self):
        n = 1;
        m = 0;
        c = 2.5;
        z = 0.5;
        type = 1;
        particle = Particle(2.5)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 3, particle)(z),
            -0.1280857374177537 + 0.01321180926927677j)

    def test_get_a_functionR3_2(self):
        n = 3;
        m = 2;
        c = 2.5;
        z = -0.5;
        type = 1;
        particle = Particle(1.5)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 3, particle)(z), -0.154763 + 0.461291j, 5)

    def test_get_a_functionR3_3(self):
        n = 5;
        m = 0;
        c = 0.5;
        z = -0.9;
        type = 1;
        particle = Particle(5)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 3, particle)(z),
            0.000642762007655166 - 0.539400871910197j)

    def test_get_a_functionR3_4(self):
        n = 4;
        m = 2;
        c = 4.5;
        z = -0.6;
        type = 1;
        psi = 6
        particle = Particle(psi)
        cv = get_cv(m, n, c, type)
        self.assertAlmostEqual(get_a_functions(m, n, c, cv, type, 3, particle)(z),
            0.03417296826676850 + 0.01608088386790631j, 5)

    def test_get_A1(self):
        particle = Particle(5)
        c1 = 0.3;
        c2 = 0.4;
        nmax = 1
        self.assertAlmostEqual(get_A(particle, c2, c1, 1, nmax)[0][0], 0.432861, 5)

    def test_get_A2(self):
        particle = Particle(3)
        c1 = 2.5;
        c2 = 1.5;
        nmax = 2
        A = get_A(particle, c2, c1, 1, nmax)
        self.assertAlmostEqual(A[0][0], 0.0450505, 5)
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[1][1], 0.238367, 5)

    def test_get_A3(self):
        particle = Particle(1.5)
        c1 = 0.5;
        c2 = 3.5;
        nmax = 3
        A = get_A(particle, c2, c1, 1, nmax)
        self.assertAlmostEqual(A[0][0], 0.0800241, 5)
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[1][1], 0.229468, 5)
        self.assertAlmostEqual(A[2][2], 0.248043, 5)
        self.assertAlmostEqual(A[0][2], -0.0150865)
        self.assertAlmostEqual(A[2][0], 0.0474487, 5)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_getA11Simplest(self):
        psi = 1.5
        c = 2
        particle = Particle(psi,d=1,eps=1)
        nmax = 3
        A = get_A11(particle,c,nmax)
        cv = get_cv(1,1,c,particle.type)
        self.assertAlmostEqual(A[0][0],rad3_cv(1,1,c,cv,particle.type,psi)[0])
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1][1],rad3_cv(1,2,c,cv,particle.type,psi)[0])
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2][2],rad3_cv(1,3,c,cv,particle.type,psi)[0])
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[0][2], 0)
        self.assertAlmostEqual(A[2][0], 0)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_getA12Simplest(self):
        psi = 3
        c = 2.5
        particle = Particle(psi,d=2.0,eps=1)
        nmax = 3
        A = get_A12(particle,c,c,nmax)
        cv = get_cv(1,1,c,particle.type)
        self.assertAlmostEqual(A[0][0],-rad1_cv(1,1,c,cv,particle.type,psi)[0])
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1][1],-rad1_cv(1,2,c,cv,particle.type,psi)[0])
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2][2],-rad1_cv(1,3,c,cv,particle.type,psi)[0])
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[0][2], 0)
        self.assertAlmostEqual(A[2][0], 0)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_getA21Simplest(self):
        psi = 3
        c = 2.5
        particle = Particle(psi,d=2.0,eps=1)
        nmax = 3
        A = get_A21(particle,c,nmax)
        cv = get_cv(1,1,c,particle.type)
        coeff = lambda nu: metric_phi(nu,particle) * metric_nu(nu,particle) / metric_psi(nu,particle) * ang1_cv(1, 1, c, cv, particle.type, nu)[0] * ang1_cv(1, 1, c, cv, particle.type, nu)[0]
        coeff = quad(coeff,-1,1)
        self.assertAlmostEqual(A[0][0],rad3_cv(1,1,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1][1],rad3_cv(1,2,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2][2],rad3_cv(1,3,c,cv,particle.type,psi)[1] * coeff)
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[0][2], 0)
        self.assertAlmostEqual(A[2][0], 0)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_getA22Simplest(self):
        psi = 3
        c = 2.5
        particle = Particle(psi,d=2.0,eps=1)
        nmax = 3
        A = get_A22(particle,c,c,nmax)
        cv = get_cv(1,1,c,particle.type)
        coeff = lambda nu: metric_phi(nu,particle) * metric_nu(nu,particle) / metric_psi(nu,particle) * ang1_cv(1, 1, c, cv, particle.type, nu)[0] * ang1_cv(1, 1, c, cv, particle.type, nu)[0]
        coeff = - quad(coeff,-1,1)
        self.assertAlmostEqual(A[0][0],rad1_cv(1,1,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1][1],rad1_cv(1,2,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2][2],rad1_cv(1,3,c,cv,particle.type,psi)[1] * coeff)
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[0][2], 0)
        self.assertAlmostEqual(A[2][0], 0)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_getA10Simplest(self):
        psi = 3
        c = 2.5
        particle = Particle(psi,d=2.0,eps=1)
        nmax = 3
        A = get_A10(particle,c,nmax)
        cv = get_cv(1,1,c,particle.type)
        self.assertAlmostEqual(A[0][0],-rad1_cv(1,1,c,cv,particle.type,psi)[0])
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1][1],-rad1_cv(1,2,c,cv,particle.type,psi)[0])
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2][2],-rad1_cv(1,3,c,cv,particle.type,psi)[0])
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[0][2], 0)
        self.assertAlmostEqual(A[2][0], 0)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_getA20Simplest(self):
        psi = 3
        c = 2.5
        particle = Particle(psi,d=2.0,eps=1)
        nmax = 3
        A = get_A20(particle,c,nmax)
        cv = get_cv(1,1,c,particle.type)
        coeff = lambda nu: metric_phi(nu,particle) * metric_nu(nu,particle) / metric_psi(nu,particle) * ang1_cv(1, 1, c, cv, particle.type, nu)[0] * ang1_cv(1, 1, c, cv, particle.type, nu)[0]
        coeff = - quad(coeff,-1,1)
        self.assertAlmostEqual(A[0][0],rad1_cv(1,1,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,2,c,particle.type)
        self.assertAlmostEqual(A[1][1],rad1_cv(1,2,c,cv,particle.type,psi)[1] * coeff)
        cv = get_cv(1,3,c,particle.type)
        self.assertAlmostEqual(A[2][2],rad1_cv(1,3,c,cv,particle.type,psi)[1] * coeff)
        self.assertAlmostEqual(A[1][0], 0)
        self.assertAlmostEqual(A[0][1], 0)
        self.assertAlmostEqual(A[0][2], 0)
        self.assertAlmostEqual(A[2][0], 0)
        self.assertAlmostEqual(A[2][1], 0)
        self.assertAlmostEqual(A[1][2], 0)

    def test_get_fullA1(self):
        c1 = 0.5;
        c2 = 3.5;
        nmax = 1;
        d = 1.5;
        psi = 1.5
        particle = Particle(psi, d)
        A = get_fullA(particle, c1, c2, nmax)
        self.assertAlmostEqual(A[0, 0], get_A(particle, c1, c1, 3, nmax))
        self.assertAlmostEqual(A[0, 1], -get_A(particle, c2, c1, 1, nmax))
        self.assertAlmostEqual(A[1, 0], get_B(particle, c1, c1, 3, nmax))
        self.assertAlmostEqual(A[1, 1], -get_C(particle, c2, c1, nmax))

    def test_get_fullA2(self):
        c1 = 0.5;
        c2 = 3.5;
        nmax = 2;
        d = 1.5;
        psi = 1.5
        particle = Particle(psi, d)
        A = get_fullA(particle, c1, c2, nmax)
        A11 = get_A(particle, c1, c1, 3, nmax)
        A12 = -get_A(particle, c2, c1, 1, nmax)
        A21 = get_B(particle, c1, c1, 3, nmax)
        A22 = -get_C(particle, c2, c1, nmax)
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

    def test_get_B1(self):
        d = 8;
        psi = 6
        c1 = 0.3;
        c2 = 0.4;
        nmax = 1
        particle = Particle(psi, d)
        self.assertAlmostEqual(get_B(particle, c2, c1, 1, nmax)[0][0], -3.69601, 5)

    def test_get_B2(self):
        d = 2;
        psi = 3
        c1 = 4;
        c2 = 5;
        nmax = 2
        particle = Particle(psi, d)
        B = get_B(particle, c2, c1, 1, nmax)
        self.assertAlmostEqual(B[0][0], 2.79893, 5)
        self.assertAlmostEqual(B[1][0], 0, 5)
        self.assertAlmostEqual(B[0][1], 0, 5)
        self.assertAlmostEqual(B[1][1], 1.55881, 5)

    def test_get_B3(self):
        d = 1;
        psi = 10
        c1 = 0.5;
        c2 = 5;
        nmax = 3
        particle = Particle(psi, d)
        B = get_B(particle, c2, c1, 1, nmax)
        self.assertAlmostEqual(B[0][0], -2.040523, 4)
        self.assertAlmostEqual(B[1][0], 0, 5)
        self.assertAlmostEqual(B[0][1], 0, 5)
        self.assertAlmostEqual(B[1][1], -4.41099, 5)
        self.assertAlmostEqual(B[1][2], 0, 5)
        self.assertAlmostEqual(B[2][1], 0, 5)
        self.assertAlmostEqual(B[0][2], 0.618216, 5)
        self.assertAlmostEqual(B[2][0], 0.373256, 5)
        self.assertAlmostEqual(B[2][2], 1.16732, 5)

    def test_getFullB1(self):
        c1 = 4.5;
        nmax = 1;
        d = 2.5;
        psi = 3.5
        particle = Particle(psi, d)
        B = get_fullB(particle, c1, nmax)
        self.assertAlmostEqual(B[0], -get_A(particle, c1, c1, 1, nmax))
        self.assertAlmostEqual(B[1], -get_B(particle, c1, c1, 1, nmax))

    def test_getFullB2(self):
        c1 = 4.5;
        nmax = 2;
        d = 2.5;
        psi = 3.5
        particle = Particle(psi, d)
        B = get_fullB(particle, c1, nmax)
        B1 = -get_A(particle, c1, c1, 1, nmax)
        B2 = -get_B(particle, c1, c1, 1, nmax)
        self.assertAlmostEqual(B[0, 0], B1[0, 0])
        self.assertAlmostEqual(B[0, 1], B1[0, 1])
        self.assertAlmostEqual(B[1, 0], B1[1, 0])
        self.assertAlmostEqual(B[1, 1], B1[1, 1])
        self.assertAlmostEqual(B[2, 0], B2[0, 0])
        self.assertAlmostEqual(B[2, 1], B2[0, 1])
        self.assertAlmostEqual(B[3, 0], B2[1, 0])
        self.assertAlmostEqual(B[3, 1], B2[1, 1])

    def test_get_C1(self):
        d = 8;
        eps = 0.7;
        psi = 6
        c1 = 0.3;
        c2 = 0.4;
        nmax = 1
        particle = Particle(psi, d, eps)
        self.assertAlmostEqual(get_C(particle, c2, c1, nmax)[0][0], -9.65321, 5)

    def test_get_C2(self):
        d = 2;
        psi = 3;
        eps = 5
        c1 = 4;
        c2 = 5;
        nmax = 2
        particle = Particle(psi, d, eps)
        C = get_C(particle, c2, c1, nmax)
        self.assertAlmostEqual(C[0][0], 0.590814, 5)
        self.assertAlmostEqual(C[1][0], 0, 5)
        self.assertAlmostEqual(C[0][1], 0, 5)
        self.assertAlmostEqual(C[1][1], 0.166018, 5)

    def test_get_C3(self):
        d = 1;
        eps = 0.3;
        psi = 10
        c1 = 0.5;
        c2 = 5;
        nmax = 3
        particle = Particle(psi, d, eps)
        C = get_C(particle, c2, c1, nmax)
        self.assertAlmostEqual(C[0][0], -6.60149, 4)
        self.assertAlmostEqual(C[1][0], 0, 5)
        self.assertAlmostEqual(C[0][1], 0, 5)
        self.assertAlmostEqual(C[1][1], -14.78301, 5)
        self.assertAlmostEqual(C[1][2], 0, 5)
        self.assertAlmostEqual(C[2][1], 0, 5)
        self.assertAlmostEqual(C[0][2], 2.00006, 5)
        self.assertAlmostEqual(C[2][0], 1.17881, 5)
        self.assertAlmostEqual(C[2][2], 3.68663, 5)
