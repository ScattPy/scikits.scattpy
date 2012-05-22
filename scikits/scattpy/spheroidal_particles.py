from scipy import *

import spheroidal

class Spheroid(object):

    def __init__(self,m,a_b,type):
        if type == 1:
            self.psi = a_b / sqrt(a_b * a_b - 1.0)
        elif type == -1:
            self.psi = 1.0 / sqrt(a_b * a_b - 1.0)
        else:
            print "Invalid value of type. Should be 1 or -1."

        self.spheroid = True

        self.m = m
        self.eps = sqrt(m)
        self.type = type
        self.k = 1.0

    def set_xl(self, xl):
        self.xl = xl
        if self.type == 1:
            self.c1 = xl / self.psi
        elif self.type == -1:
            self.c1 = xl / sqrt(self.psi * self.psi + 1.0)
        else:
            print "Invalid value of type. Should be 1 or -1. Current value is "+str(self.type)
        self.xv = self.c1 * pow(self.psi * (self.psi * self.psi - self.type),1.0/3.0)
        self.initialize()

    def set_xv(self,xv):
        self.xv = xv
        self.c1 = xv / pow(self.psi * (self.psi * self.psi - self.type),1.0/3.0)
        if self.type == 1:
            self.xl = self.c1 * self.psi
        elif self.type == -1:
            self.xl = self.c1 * sqrt(self.psi * self.psi + 1.0)
        else:
            print "Invalid value of type. Should be 1 or -1."
        self.initialize()

    def initialize(self):
        self.c2 = self.m * self.c1
        self.d = 2.0 * self.c1

    def function(self, nu):
        return self.psi

    def derivative(self, nu):
        return 0.0

#Incedent waves
class TEInputWave:
    alpha = 0

    def __init__(self, alpha):
        self.alpha = alpha

    #according to formulas 52
    def getB(self, particle, l):
        return 0

    def getA(self, particle, c, l):
        cv = spheroidal.get_cv(1, l, c, particle.type)
        return -2 * pow(1j, l) * spheroidal.ang1_cv(1, l, c, cv, particle.type, cos([self.alpha]))[0]


class TMInputWave:
    alpha = 0

    def __init__(self, alpha):
        self.alpha = alpha

    #according to formulas 53
    def getB(self, particle, c, l):
        cv = spheroidal.get_cv(1, l, c, particle.type)
        return 2 * pow(1j, l) * spheroidal.ang1_cv(1, l, c, cv, particle.type, cos([self.alpha]))[0]

    def getA(self, particle, l):
        return 0.0
