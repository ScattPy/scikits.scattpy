import scipy
import scipy.integrate
import scipy.linalg

from scipy import *

from spheroidal_functions import *
from spheroidal_svm import *

#quad integration for complex numbers
from spheroidal_svm import get_fullA

def quad(func, a, b, **kwargs):
    def real_func(x):
        return scipy.real(func(x))
    def imag_func(x):
        return scipy.imag(func(x))
    real_integral = scipy.integrate.quad(real_func, a, b, **kwargs)
    imag_integral = scipy.integrate.quad(imag_func, a, b, **kwargs)
    return real_integral[0] + 1j*imag_integral[0]

# ---- Generation of A matrices

#according to (86)
def get_A11(particle, c1, nmax):
    return get_A(particle, c1, c1, 3, nmax).transpose()


def get_A12(particle, c1, c2, nmax):
    return -get_A(particle, c2, c1, 1, nmax).transpose()


def get_A10(particle, c1, nmax):
    return -get_A(particle, c1, c1, 1, nmax).transpose()


def get_A21(particle, c1, nmax):
    return get_B(particle, c1, c1, 3, nmax).transpose()


def get_A22(particle, c1, c2, nmax):
    return -get_C(particle, c2, c1, nmax).transpose()


def get_A20(particle, c1, nmax):
    return -get_B(particle, c1, c1, 1, nmax).transpose()

#according to (81)
def get_a_functions(m, n, c, cv, type, rank, particle):
    if rank == 1:
        return lambda nu: rad1_cv(m, n, c, cv, type, particle.function(nu))[0] * ang1_cv(m, n, c, cv, type, nu)[0]
    elif rank == 3:
        return lambda nu: rad3_cv(m, n, c, cv, type, particle.function(nu))[0] * ang1_cv(m, n, c, cv, type, nu)[0]

#according to (81)
def get_b_functions(m, n, c, cv, type, rank, particle):
    if rank == 1:
        return lambda nu: (metric_nu(nu, particle) / metric_psi(nu, particle)
                           * rad1_cv(m, n, c, cv, type, particle.function(nu))[1] * ang1_cv(m, n, c, cv, type, nu)[0]
                           - particle.derivative(nu) * metric_psi(nu, particle) / metric_nu(nu, particle)
                             * rad1_cv(m, n, c, cv, type, particle.function(nu))[0] * ang1_cv(m, n, c, cv, type, nu)[1])
    elif rank == 3:
        return lambda nu: (metric_nu(nu, particle) / metric_psi(nu, particle)
                           * rad3_cv(m, n, c, cv, type, particle.function(nu))[1] * ang1_cv(m, n, c, cv, type, nu)[0]
                           - particle.derivative(nu) * metric_psi(nu, particle) / metric_nu(nu, particle)
                             * rad3_cv(m, n, c, cv, type, particle.function(nu))[0] * ang1_cv(m, n, c, cv, type, nu)[1])

#according to (82)
def get_c_functions(m, n, c, cv, type, rank, particle):
    eps = particle.eps
    #according to (30)
    delta = lambda nu: metric_phi(nu, particle)
    #according to (28)
    return lambda nu: get_b_functions(m, n, c, cv, type, 1, particle)(nu) / eps -\
                      (1.0 / eps - 1.0) * IzIt(nu,particle) / delta(nu) * get_a_functions(m, n, c, cv, type, 1, particle)(nu) \
                        * get_integral_metric(particle)(nu)

#-------------------------Metric coefficients ------------------------------------------------

#according to (81)
def get_integral_metric(particle):
    return lambda nu: sqrt(metric_nu(nu, particle) * metric_nu(nu, particle)
        + particle.derivative(nu) * particle.derivative(nu) * metric_psi(nu, particle) * metric_psi(nu, particle))

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

def IzIt(nu,particle):
    return particle.d / 2.0 * (particle.derivative(nu) * nu + particle.function(nu))\
    / get_integral_metric(particle)(nu)

def IzIn(nu,particle):
    return particle.d / 2.0 * (nu * metric_nu(nu,particle) / metric_psi(nu,particle) -
                               particle.derivative(nu) * particle.function(nu) * metric_psi(nu,particle) / metric_nu(nu,particle))\
    / get_integral_metric(particle)(nu)

def RIn(nu,particle):
    return (particle.d / 2.0)**2 * (particle.function(nu) * metric_nu(nu,particle) / metric_psi(nu,particle) -
                                    particle.type * nu *particle.derivative(nu) * metric_psi(nu,particle) / metric_nu(nu,particle)) \
                                    / get_integral_metric(particle)(nu)
def RIt(nu,particle):
    return (particle.d / 2.0)**2 * (particle.derivative(nu)*particle.function(nu)+particle.type*nu) / get_integral_metric(particle)(nu)
#----------------------------------------------------------------------------------------------

#according to (83)
def get_A(particle, c2, c1, rank, nmax):
    A = zeros((nmax, nmax),dtype=complex)
    type = particle.type
    m = 1
    for i in range(nmax):
        for k in range(nmax):
            l = i + m
            n = k + m
            cv_l = get_cv(m, l, c2, type)
            cv_n = get_cv(m, n, c1, type)
            func = lambda nu: get_a_functions(m, l, c2, cv_l, type, rank, particle)(nu) * ang1_cv(m, n, c1, cv_n, particle.type, nu)[0]
            A[i][k] = quad(func, -1, 1)
    return A

#according to (84)
def get_Z(get_z_functions, particle, c2, c1, rank, nmax):
    Z = zeros((nmax, nmax),dtype=complex)
    type = particle.type
    m = 1
    for i in range(nmax):
        for k in range(nmax):
            l = i + m
            n = k + m
            cv_l = get_cv(m, l, c2, type)
            cv_n = get_cv(m, n, c1, type)
            func = lambda nu: get_z_functions(m, l, c2, cv_l, type, rank, particle)(nu) *\
                              ang1_cv(m, n, c1, cv_n, particle.type, nu)[0] * metric_phi(nu, particle)
            Z[i][k] = quad(func, -1, 1)
    return Z


def get_C(particle, c2, c1, nmax):
    return get_Z(get_c_functions, particle, c2, c1, 1, nmax)


def get_B(particle, c2, c1, rank, nmax):
    return get_Z(get_b_functions, particle, c2, c1, rank, nmax)

#-------Solve equation and find solution of scattering

def get_Bin(inputWave, particle, c, nmax):
    b = zeros((nmax, 1),dtype=complex)
    for i in range(nmax):
        l = i + 1
        b[i] = inputWave.getB(particle, c, l)
    return b

#Return b_sca and b_int. b_sca = result[0] and b_int = result[1]
def getSolution(particle, inputWave, c1, c2, nmax):
    A = get_fullA(particle, c1, c2, nmax)
    B = get_fullB(particle, c1, nmax) * get_Bin(inputWave, particle, c1, nmax)
    x = scipy.linalg.solve(A, B)
    return (x[0:nmax + 1], x[nmax + 1:])