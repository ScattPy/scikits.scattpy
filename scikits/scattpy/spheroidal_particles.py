from scipy import sqrt

class ProlateSpheroid:
    type = 1
    d = 1
    eps = 1

    def function(self, nu):
        return 1.1

    def derivative(self, nu):
        return 0


def metric_psi(nu, particle):
    psi = particle.function(nu)
    return particle.d / 2.0 * sqrt((psi * psi - particle.type * nu * nu) / (psi * psi - particle.type))


def metric_nu(nu, particle):
    psi = particle.function(nu)
    return particle.d / 2.0 * sqrt((psi * psi - particle.type * nu * nu) / (1 - nu * nu))


def metric_phi(nu, particle):
    psi = particle.function(nu)
    return particle.d / 2.0 * sqrt((psi * psi - particle.type) * (1 - nu * nu))


class PlainInputWave:
    def getB(self, m, l):
        return 0

    def getA(self, m, l):
        return 0


