from numpy import zeros, asarray, mat
import scipy

from spheroidal_functions import *

import spheroidal
from spheroidal_svm import SpheroidalSVM


class SpheroidalEBCM:
    def __init__(self, particle, c2, c1, nmax):
        self.particle = particle
        self.c2 = c2
        self.c1 = c1
        self.nmax = nmax

    #according to (88)
    def get_BSi(self):
        Bs = zeros((self.nmax, self.nmax), dtype=complex)
        type = self.particle.type
        m = 1
        for i in range(self.nmax):
            for k in range(self.nmax):
                n = i + m
                l = k + m
                cv_l = get_cv(m, l, self.c2, type)
                cv_n = get_cv(m, n, self.c1, type)
                func = lambda nu: (spheroidal.get_a_function(m, l, self.c2, cv_l, 1, self.particle)(nu) *\
                                  spheroidal.get_b_function(m, n, self.c1, cv_n, 3, self.particle)(nu) -\
                                  spheroidal.get_c_function(m, l, self.c2, cv_l, 1, self.particle)(nu) *\
                                  spheroidal.get_a_function(m, n, self.c1, cv_n, 3, self.particle)(nu)) * spheroidal.metric_phi(nu, self.particle) *\
                                  spheroidal.get_integral_metric(self.particle)(nu)

                Bs[i][k] = spheroidal.quad(func, -1, 1)
        return  1j * mat(Bs)

    #according to (88)
    def get_BRi(self):
        Br = zeros((self.nmax, self.nmax), dtype=complex)
        type = self.particle.type
        m = 1
        for i in range(self.nmax):
            for k in range(self.nmax):
                n = i + m
                l = k + m
                cv_l = get_cv(m, l, self.c2, type)
                cv_n = get_cv(m, n, self.c1, type)
                func = lambda nu: (spheroidal.get_a_function(m, l, self.c2, cv_l, 1, self.particle)(nu) * \
                                  spheroidal.get_b_function(m, n, self.c1, cv_n, 1, self.particle)(nu) - \
                                  spheroidal.get_c_function(m, l, self.c2, cv_l, 1, self.particle)(nu) * \
                                  spheroidal.get_a_function(m, n, self.c1, cv_n, 1, self.particle)(nu)) * spheroidal.metric_phi(nu, self.particle) *\
                                  spheroidal.get_integral_metric(self.particle)(nu)

                Br[i][k] = spheroidal.quad(func, -1, 1)
        return  -1j * mat(Br)

    #according to (99)
    def get_BSm(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        Bs = 1j * ( svm.get_B(3) * svm.get_A(self.c2, self.c1, 1).T -
                    svm.get_A(self.c1, self.c1, 3) * svm.get_C().T)
        return Bs

    #according to (99)
    def get_BRm(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        Br = -1j * ( svm.get_B(1) * svm.get_A(self.c2, self.c1, 1).T -
                    svm.get_A(self.c1, self.c1, 1) * svm.get_C().T)
        return Br


    #Return b_sca and b_int. b_sca = result[0] and b_int = result[1]
    def getSolution(self, inputWave):
        b_in = spheroidal.get_Bin(inputWave, self.particle, self.c1, self.nmax)
        b_int = scipy.linalg.solve(self.get_BSi(), -b_in)
        b_sca = asarray(self.get_BRi() * b_int)
        return (b_sca, b_int)

    def getMatrixSolution(self, inputWave):
        b_in = spheroidal.get_Bin(inputWave, self.particle, self.c1, self.nmax)
        b_int = scipy.linalg.solve(self.get_BSm(), -b_in)
        b_sca = asarray(self.get_BRm() * b_int)
        return (b_sca, b_int)
