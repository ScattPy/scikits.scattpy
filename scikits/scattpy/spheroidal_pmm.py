from numpy import zeros, bmat, mat, conjugate
import scipy

from spheroidal_functions import *

import spheroidal
from spheroidal_svm import SpheroidalSVM


class SpheroidalPMM:

    def __init__(self,particle,c2,c1,nmax):
        self.particle=particle
        self.c2 = c2
        self.c1 = c1
        self.nmax = nmax

    def get_A(self,function):
        A = zeros((self.nmax, self.nmax), dtype=complex)
        m = 1
        for i in range(self.nmax):
            for k in range(self.nmax):
                n = i + m
                l = k + m
                func = lambda nu: function(m, n,l)(nu) * \
                                   spheroidal.metric_phi(nu, self.particle) *\
                                   spheroidal.get_integral_metric(self.particle)(nu)

                A[i][k] = spheroidal.quad(func, -1, 1)
        return  mat(A)

    # ---- Generation of A matrices

    def function_a11(self,m, n, l):
        cv_l = get_cv(m, l, self.c1, self.particle.type)
        cv_n = get_cv(m, n, self.c1, self.particle.type)
        return lambda nu: conjugate(spheroidal.get_a_function(m,n,self.c1,cv_n,3,self.particle)(nu)) *\
                          spheroidal.get_a_function(m,l,self.c1,cv_l,3,self.particle)(nu) +\
                          conjugate(spheroidal.get_b_function(m,n,self.c1,cv_n,3,self.particle)(nu)) *\
                          spheroidal.get_b_function(m,l,self.c1,cv_l,3,self.particle)(nu)

    def function_a12(self,m, n, l):
        cv_l = get_cv(m, l, self.c2, self.particle.type)
        cv_n = get_cv(m, n, self.c1, self.particle.type)
        return lambda nu: -(conjugate(spheroidal.get_a_function(m,n,self.c1,cv_n,3,self.particle)(nu)) *\
                          spheroidal.get_a_function(m,l,self.c2,cv_l,1,self.particle)(nu) +\
                          conjugate(spheroidal.get_b_function(m,n,self.c1,cv_n,3,self.particle)(nu)) *\
                          spheroidal.get_c_function(m,l,self.c2,cv_l,1,self.particle)(nu))

    def function_a21(self,m, n, l):
        cv_l = get_cv(m, l, self.c1, self.particle.type)
        cv_n = get_cv(m, n, self.c2, self.particle.type)
        return lambda nu: -(conjugate(spheroidal.get_a_function(m,n,self.c2,cv_n,1,self.particle)(nu)) *\
                          spheroidal.get_a_function(m,l,self.c1,cv_l,3,self.particle)(nu) +\
                          conjugate(spheroidal.get_c_function(m,n,self.c2,cv_n,1,self.particle)(nu)) *\
                          spheroidal.get_b_function(m,l,self.c1,cv_l,3,self.particle)(nu))

    def function_a22(self,m, n, l):
        cv_l = get_cv(m, l, self.c2, self.particle.type)
        cv_n = get_cv(m, n, self.c2, self.particle.type)
        return lambda nu: conjugate(spheroidal.get_a_function(m,n,self.c2,cv_n,1,self.particle)(nu)) *\
                          spheroidal.get_a_function(m,l,self.c2,cv_l,1,self.particle)(nu) +\
                          conjugate(spheroidal.get_c_function(m,n,self.c2,cv_n,1,self.particle)(nu)) *\
                          spheroidal.get_c_function(m,l,self.c2,cv_l,1,self.particle)(nu)

    def function_a10(self,m, n, l):
        cv_l = get_cv(m, l, self.c1, self.particle.type)
        cv_n = get_cv(m, n, self.c1, self.particle.type)
        return lambda nu: -(conjugate(spheroidal.get_a_function(m,n,self.c1,cv_n,3,self.particle)(nu)) *\
                          spheroidal.get_a_function(m,l,self.c1,cv_l,1,self.particle)(nu) +\
                          conjugate(spheroidal.get_b_function(m,n,self.c1,cv_n,3,self.particle)(nu)) *\
                          spheroidal.get_b_function(m,l,self.c1,cv_l,1,self.particle)(nu))

    def function_a20(self,m, n, l):
        cv_l = get_cv(m, l, self.c1, self.particle.type)
        cv_n = get_cv(m, n, self.c2, self.particle.type)
        return lambda nu: conjugate(spheroidal.get_a_function(m,n,self.c2,cv_n,1,self.particle)(nu)) *\
                          spheroidal.get_a_function(m,l,self.c1,cv_l,1,self.particle)(nu) +\
                          conjugate(spheroidal.get_c_function(m,n,self.c2,cv_n,1,self.particle)(nu)) *\
                          spheroidal.get_b_function(m,l,self.c1,cv_l,1,self.particle)(nu)

    def get_A11i(self):
        return self.get_A(lambda m,n,l: self.function_a11(m,n,l))

    def get_A12i(self):
        return self.get_A(lambda m,n,l: self.function_a12(m,n,l))

    def get_A21i(self):
        return self.get_A(lambda m,n,l: self.function_a21(m,n,l))

    def get_A22i(self):
        return self.get_A(lambda m,n,l: self.function_a22(m,n,l))

    def get_A10i(self):
        return self.get_A(lambda m,n,l: self.function_a10(m,n,l))

    def get_A20i(self):
        return self.get_A(lambda m,n,l: self.function_a20(m,n,l))

    #according to (101)
    def get_A11m(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        return conjugate(svm.get_A(self.c1,self.c1,3)) * svm.get_A(self.c1,self.c1,3).T + \
            conjugate(svm.get_B(3)) * svm.get_B(3).T


    def get_A12m(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        return -(conjugate(svm.get_A(self.c1,self.c1,3)) * svm.get_A(self.c2,self.c1,1).T +\
               conjugate(svm.get_B(3)) * svm.get_C().T)


    def get_A21m(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        return -(conjugate(svm.get_A(self.c2,self.c1,1)) * svm.get_A(self.c1,self.c1,3).T +\
               conjugate(svm.get_C()) * svm.get_B(3).T)

    def get_A22m(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        return conjugate(svm.get_A(self.c2,self.c1,1)) * svm.get_A(self.c2,self.c1,1).T +\
               conjugate(svm.get_C()) * svm.get_C().T

    def get_A10m(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        return -(conjugate(svm.get_A(self.c1,self.c1,3)) * svm.get_A(self.c1,self.c1,1).T +\
               conjugate(svm.get_B(3)) * svm.get_B(1).T)

    def get_A20m(self):
        svm = SpheroidalSVM(self.particle, self.c2, self.c1, self.nmax)
        return conjugate(svm.get_A(self.c2,self.c1,1)) * svm.get_A(self.c1,self.c1,1).T +\
                 conjugate(svm.get_C()) * svm.get_B(1).T

    def get_fullBm(self):
        return bmat([[self.get_A10m()], [self.get_A20m()]])

    def get_fullAm(self):
        A11 = self.get_A11m()
        A12 = self.get_A12m()
        A21 = self.get_A21m()
        A22 = self.get_A22m()
        return bmat([[A11, A12], [A21, A22]])

    def getMatrixSolution(self,inputWave):
        A = self.get_fullAm()
        B = self.get_fullBm() * spheroidal.get_Bin(inputWave, self.particle, self.c1, self.nmax)
        x = -scipy.linalg.solve(A, B)
        return (x[0:self.nmax + 1], x[self.nmax + 1:])

    def get_fullBi(self):
        return bmat([[self.get_A10i()], [self.get_A20i()]])

    def get_fullAi(self):
        A11 = self.get_A11i()
        A12 = self.get_A12i()
        A21 = self.get_A21i()
        A22 = self.get_A22i()
        return bmat([[A11, A12], [A21, A22]])

    def getSolution(self,inputWave):
        A = self.get_fullAi()
        B = self.get_fullBi() * spheroidal.get_Bin(inputWave, self.particle, self.c1, self.nmax)
        x = -scipy.linalg.solve(A, B)
        return (x[0:self.nmax + 1], x[self.nmax + 1:])