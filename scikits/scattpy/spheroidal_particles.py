from scipy import sqrt
from numpy import cos

import spheroidal

class ProlateSpheroid:
    type = 1
    d = 1
    eps = 1
    c = 1

    def function(self, nu):
        return 1.5

    def derivative(self, nu):
        return 0

#according to formula 25
def metric_psi(nu, particle):
    psi = particle.function(nu)
    return particle.d / 2.0 * sqrt((psi * psi - particle.type * nu * nu) / (psi * psi - particle.type))


def metric_nu(nu, particle):
    psi = particle.function(nu)
    return particle.d / 2.0 * sqrt((psi * psi - particle.type * nu * nu) / (1 - nu * nu))


def metric_phi(nu, particle):
    psi = particle.function(nu)
    return particle.d / 2.0 * sqrt((psi * psi - particle.type) * (1 - nu * nu))

#Incedent waves
class TEInputWave:
    alpha = 0

    def __init__(self,alpha):
        self.alpha = alpha

    #according to formulas 52
    def getB(self, particle, l):
        return 0

    def getA(self, particle, l):
        c = particle.c
        cv = spheroidal.get_cv(1, l, c, particle.type)
        return -2 * pow(1j,l) * spheroidal.ang1_cv(1,l,c,cv,particle.type,cos([self.alpha]))[0] / spheroidal.get_norm_cv(1,l,c, cv, particle.type)


class TMInputWave:
    alpha = 0

    def __init__(self,alpha):
        self.alpha = alpha

    #according to formulas 53
    def getB(self,particle, l):
        c = particle.c
        cv = spheroidal.get_cv(1, l, c, particle.type)
        return 2 * pow(1j,l) * spheroidal.ang1_cv(1,l,c,cv,particle.type,cos([self.alpha]))[0] / spheroidal.get_norm_cv(1,l,c, cv, particle.type)

    def getA(self, particle, l):
        return 0.0
