import numpy
import unittest
import scipy.linalg

from scipy.special import *

from spheroidal import *
from spheroidal_particles import Spheroid

class testSpheroidalProlateNorm(unittest.TestCase):
    numpy.seterr('print')
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
    def test_axinorm1(self):
        c=1; m=1; n=5
        self.check_norm(m,n,c)

    def test_axinorm5(self):
        c=5; m=1; n=10
        self.check_norm(m,n,c)

    def test_axinorm5(self):
        c=10; m=1; n=20
        self.check_norm(m,n,c)

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
        func = lambda x: power(pro_ang1(m, n, c, x)[0], 2)
        self.assertAlmostEqual((quad(func, -1, 1, limit = 300, epsabs=10e-10,epsrel=10e-10))/(norm*norm), 1, places = 7)

class testSpheroidalOblateNorm(testSpheroidalProlateNorm):
    def check_norm(self, m, n, c):
        norm = get_obl_norm(m,n,c)
        func = lambda x: power(obl_ang1(m, n, c, x)[0] / norm, 2)
        self.assertAlmostEqual((quad(func, -1, 1, limit = 300,epsabs=10e-10,epsrel=10e-10)), 1, places = 7)

class testAngularRelationsProlate(unittest.TestCase):
    def testEven1(self):
        m = 1; n = 2; c = 2; x = 0.5
        self.assertAlmostEqual(pro_ang1(m, n, c, x)[0],pow(-1,n-m) * pro_ang1(m, n, c, -x)[0])
        self.assertAlmostEqual(obl_ang1(m, n, c, x)[0],pow(-1,n-m) * obl_ang1(m, n, c, -x)[0])

    def testEven2(self):
        m = 1; n = 3; c = 2; x = -0.5
        self.assertAlmostEqual(pro_ang1(m, n, c, x)[0],pow(-1,n-m) * pro_ang1(m, n, c, -x)[0])
        self.assertAlmostEqual(obl_ang1(m, n, c, x)[0],pow(-1,n-m) * obl_ang1(m, n, c, -x)[0])

    def testSimplestCaseProlate1(self):
        m = 1; n = 2; c = numpy.pi * n / 2; type = 1; x = 0.5
        cv = get_cv(m, n, c, type)
        norm = get_norm_cv(m, n, c, cv, type)
        self.assertAlmostEqual(ang1_cv(m, n, c,cv,type, x)[0] * norm, - 2 / (n * numpy.pi) * scipy.special.lpmn(1,n,0)[1][m][n] /sqrt(1-x*x) * sin(numpy.pi * n * x/2))

    def testSimplestCaseProlate2(self):
        m = 1; n = 3; c = numpy.pi * n / 2; type = 1; x = 0.3
        cv = get_cv(m, n, c, type)
        norm = get_norm_cv(m, n, c, cv, type)
        self.assertAlmostEqual(ang1_cv(m, n, c,cv,type, x)[0] * norm,  - scipy.special.lpmv(1,n,0) /sqrt(1-x*x) * cos(numpy.pi * n * x/2))

    def testSimplestCaseProlate3(self):
        m = 1; n = 4; c = numpy.pi * n / 2; type = 1; x = -0.5
        cv = get_cv(m, n, c, type)
        norm = get_norm_cv(m, n, c, cv, type)
        self.assertAlmostEqual(ang1_cv(m, n, c,cv,type, x)[0] * norm, 2 / (n * numpy.pi) * scipy.special.lpmn(1,n,0)[1][m][n] /sqrt(1-x*x) * sin(numpy.pi * n * x/2))

    def testSimplestCaseProlate4(self):
        m = 1; n = 5; c = numpy.pi * n / 2; type = 1; x = -0.3
        cv = get_cv(m, n, c, type)
        norm = get_norm_cv(m, n, c, cv, type)
        self.assertAlmostEqual(ang1_cv(m, n, c,cv,type, x)[0] * norm,  - scipy.special.lpmv(1,n,0) /sqrt(1-x*x) * cos(numpy.pi * n * x/2))

    def testSimplestCaseProlate5(self):
        m = 1; n = 6; c = numpy.pi * n / 2; type = 1; x = -0.3
        cv = get_cv(m, n, c, type)
        norm = get_norm_cv(m, n, c, cv, type)
        self.assertAlmostEqual(ang1_cv(m, n, c,cv,type, x)[0] * norm, - 2 / (n * numpy.pi) * scipy.special.lpmn(1,n,0)[1][m][n] /sqrt(1-x*x) * sin(numpy.pi * n * x/2))

    def testSimplestCaseProlate6(self):
        m = 1; n = 7; c = numpy.pi * n / 2; type = 1; x = -0.3
        cv = get_cv(m, n, c, type)
        norm = get_norm_cv(m, n, c, cv, type)
        self.assertAlmostEqual(ang1_cv(m, n, c,cv,type, x)[0] * norm,  - scipy.special.lpmv(1,n,0) /sqrt(1-x*x) * cos(numpy.pi * n * x/2))

class testRadialRelations(unittest.TestCase):

    def testR1R2Wronskian1(self):
        m=1;n=3;c = 2; psi =3.0
        [R1,dR1] = scipy.special.pro_rad1(m,n,c,psi)
        [R2,dR2] = scipy.special.pro_rad2(m,n,c,psi)
        W = scipy.mat([[R1, dR1],[R2, dR2]])
        self.assertAlmostEqual(scipy.linalg.det(W),1/(c*(psi**2-1)))

    def testR1R2Wronskian2(self):
        m=1;n=2;c = 3.0; psi = 2.5
        [R1,dR1] = scipy.special.obl_rad1(m,n,c,psi)
        [R2,dR2] = scipy.special.obl_rad2(m,n,c,psi)
        W = scipy.mat([[R1, dR1],[R2, dR2]])
        self.assertAlmostEqual(scipy.linalg.det(W),1/(c*(psi**2+1)))

    def testR1R3Wronskian1(self):
        m=1;n=3;c = 2; psi =3.0
        type = 1
        cv = get_cv(m, n, c, type)
        [R1,dR1] = scipy.special.pro_rad1(m,n,c,psi)
        [R3,dR3] = rad3_cv(m,n,c,cv,type,psi)
        W = scipy.mat([[R1, dR1],[R3, dR3]])
        self.assertAlmostEqual(scipy.linalg.det(W),1j/(c*(psi**2-1)))

    def testR1R3Wronskian2(self):
        m=1;n=2;c = 3.0; psi = 2.5
        type = -1
        cv = get_cv(m, n, c, type)
        [R1,dR1] = scipy.special.obl_rad1(m,n,c,psi)
        [R3,dR3] = rad3_cv(m,n,c,cv,type,psi)
        W = scipy.mat([[R1, dR1],[R3, dR3]])
        self.assertAlmostEqual(scipy.linalg.det(W),1j/(c*(psi**2+1)))



class testMetricCoefficients(unittest.TestCase):

    def testDelta1(self):
        m=1.5;a_b=3;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nu = 0.5
        self.assertAlmostEquals(IzIn(nu,particle)*RIt(nu,particle)-IzIt(nu,particle)*RIn(nu,particle),-metric_phi(nu,particle))

    def testDelta2(self):
        m=1.5;a_b=3;type=-1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nu = -0.3
        self.assertAlmostEquals(IzIn(nu,particle)*RIt(nu,particle)-IzIt(nu,particle)*RIn(nu,particle),-metric_phi(nu,particle))

    def testDelta3(self):
        m=2;a_b=2;type=1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nu = -0.5
        self.assertAlmostEquals(IzIn(nu,particle)*RIt(nu,particle)-IzIt(nu,particle)*RIn(nu,particle),-metric_phi(nu,particle))

    def testDelta4(self):
        m=2;a_b=2;type=-1
        particle = Spheroid(m,a_b,type)
        particle.set_xl(1.5)
        nu = 0.3
        self.assertAlmostEquals(IzIn(nu,particle)*RIt(nu,particle)-IzIt(nu,particle)*RIn(nu,particle),-metric_phi(nu,particle))

class testDeltaAngularFunctions(unittest.TestCase):
    def testDelta_n_equal(self):
        m=1;
        c1 = 2; n1 = 3
        c2 = 2; n2 = 3;
        norm1 = get_pro_norm(m,n1,c1)
        norm2 = get_pro_norm(m,n2,c2)
        func = lambda x: pro_ang1(m, n1, c1, x)[0] / norm1 * pro_ang1(m, n2, c2, x)[0] / norm2
        self.assertAlmostEqual((quad(func, -1, 1)), 1, places = 7)

    def testDelta_n_different1(self):
        m=1;
        c1 = 2; n1 = 4
        c2 = 2; n2 = 3;
        func = lambda x: pro_ang1(m, n1, c1, x)[0]  * pro_ang1(m, n2, c2, x)[0]
        self.assertAlmostEqual((quad(func, -1, 1)), 0, places = 7)

    def testDelta_n_different2(self):
        m=1; type = 1
        c1 = 10; n1 = 1
        c2 = 10; n2 = 3
        cv1 = get_cv(m,n1,c1,type)
        cv2 = get_cv(m,n2,c2,type)
        func = lambda x: ang1_cv(m, n1, c1,cv1,type, x)[0]  * ang1_cv(m, n2, c2,cv2,type, x)[0]
        self.assertAlmostEqual((quad(func, -1, 1,epsabs=1e-12,epsrel=1e-12)), 0, places = 10)

    def testDelta_n_different(self):
        m=1; type = 1
        c1 = 2; n1 = 4
        c2 = 2; n2 = 2;
        cv1 = get_cv(m,n1,c1,type)
        cv2 = get_cv(m,n2,c2,type)
        func = lambda x: ang1_cv(m, n1, c1,cv1,type, x)[0]  * ang1_cv(m, n2, c2,cv2,type, x)[0]
        self.assertAlmostEqual((quad(func, -1, 1)), 0, places = 7)

    def testDelta_simmetric(self):
        m=1;
        c1 = 2; n1 = 4
        c2 = 3; n2 = 4;
        func = lambda x: pro_ang1(m, n1, c1, x)[0]  * numpy.conjugate(pro_ang1(m, n2, c2, x)[0])
        func2 = lambda x: pro_ang1(m, n2, c2, x)[0]  * numpy.conjugate(pro_ang1(m, n1, c1, x)[0])
        self.assertAlmostEqual(quad(func, -1, 1), quad(func2, -1, 1), places = 7)

class testFactorial(unittest.TestCase):
    
    def test_factorial(self):
        self.assertEqual(get_norm_factorial(1,0), 1)
        self.assertEqual(get_norm_factorial(5,0), 1)
        self.assertEqual(get_norm_factorial(0,0), 1)

        self.assertEqual(get_norm_factorial(3,1), 20)
        self.assertEqual(get_norm_factorial(2,2), 360)


if __name__ == '__main__':
    unittest.main()
