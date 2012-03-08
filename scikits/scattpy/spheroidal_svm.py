from numpy import zeros, bmat, mat
import scipy
import scipy.linalg

from spheroidal_functions import *

import spheroidal


class SpheroidalSVM:

    def __init__(self,particle,c2,c1,nmax):
        self.particle=particle
        self.c2 = c2
        self.c1 = c1
        self.nmax = nmax

    #according to (83)
    def get_A(self, c2,c1,rank):
        A = zeros((self.nmax, self.nmax),dtype=complex)
        type = self.particle.type
        m = 1
        for i in range(self.nmax):
            for k in range(self.nmax):
                l = i + m
                n = k + m
                cv_l = get_cv(m, l, c2, type)
                cv_n = get_cv(m, n, c1, type)
                func = lambda nu: spheroidal.get_a_function(m, l, c2, cv_l, rank, self.particle)(nu) * ang1_cv(m, n, c1, cv_n, type, nu)[0]
                A[i][k] = spheroidal.quad(func, -1, 1)
        return mat(A)

    #according to (84)
    def __get_Z(self,z_function,c2,c1, rank):
        Z = zeros((self.nmax, self.nmax),dtype=complex)
        type = self.particle.type
        m = 1
        for i in range(self.nmax):
            for k in range(self.nmax):
                l = i + m
                n = k + m
                cv_l = get_cv(m, l, c2, type)
                cv_n = get_cv(m, n, c1, type)
                func = lambda nu: z_function(m, l, c2, cv_l, rank, self.particle)(nu) *\
                                  ang1_cv(m, n, c1, cv_n, type, nu)[0] * spheroidal.metric_phi(nu, self.particle) * \
                                  spheroidal.get_integral_metric(self.particle)(nu)
                Z[i][k] = spheroidal.quad(func, -1, 1)
        return mat(Z)


    def get_C(self):
        return self.__get_Z(spheroidal.get_c_function, self.c2, self.c1, 1)


    def get_B(self, rank):
        return self.__get_Z(spheroidal.get_b_function, self.c1, self.c1, rank)

    # ---- Generation of A matrices

    #according to (86)
    def get_A11(self):
        return self.get_A(self.c1, self.c1, 3).T


    def get_A12(self):
        return -self.get_A(self.c2, self.c1, 1).T


    def get_A10(self):
        return -self.get_A(self.c1, self.c1, 1).T


    def get_A21(self):
        return self.get_B(3).T


    def get_A22(self):
        return -self.get_C().T


    def get_A20(self):
        return -self.get_B(1).T

    #according to (85)
    def get_fullB(self):
        return bmat([[self.get_A10()], [self.get_A20()]])

    def get_fullA(self):
        A11 = self.get_A11()
        A12 = self.get_A12()
        A21 = self.get_A21()
        A22 = self.get_A22()
        return bmat([[A11, A12], [A21, A22]])

    def getSolution(self,inputWave):
        A = self.get_fullA()
        B = self.get_fullB() * spheroidal.get_Bin(inputWave, self.particle, self.c1, self.nmax)
        x = -scipy.linalg.solve(A, B)
        return (x[0:self.nmax + 1], x[self.nmax + 1:])