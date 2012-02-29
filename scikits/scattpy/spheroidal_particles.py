from numpy import cos

import spheroidal

class ProlateSpheroid:
    type = 1
    d = 4.
    c = 5.

    def __init__(self,psi=2.,eps=1.,c=1,derivative=0):
        self.psi = psi
        self.deriv = derivative
        self.eps = eps
        self.c = c

    def function(self, nu):
        return self.psi

    def derivative(self, nu):
        return self.deriv


class OblateSpheroid:
    type = -1
    d = 1
    eps = 1
    c = 1

    def __init__(self,psi=2,derivative=0):
        self.psi = psi
        self.deriv = derivative

    def function(self, nu):
        return self.psi

    def derivative(self, nu):
        return self.deriv

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
