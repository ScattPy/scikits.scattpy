import scipy.special as special
import scipy.linalg as linalg
from scipy.integrate import quad
import f_spheroid

from scipy import *

# ---- Calculation of norm for spheroidal angular functions ----
from spheroidal_particles import metric_phi, metric_nu, metric_psi

def get_norm_cv(m, n, c, cv, type):
    eps = 1E-20
    value = 1
    sum = 0
    k = 0
    d = f_spheroid.sdmn(m, n, c, cv, type)

    isEven = (n - m) % 2
    if(isEven == 1):
        k = 1

    while (value > eps):
        if (type == 1):
            value = d[k / 2] * d[k / 2] / (2.0 * k + 2.0 * m + 1.0) * get_norm_factorial(k, m)
        else:
            value = d[k / 2] * d[k / 2] / (2.0 * k + 2.0 * m + 1.0) * get_norm_factorial(k, m)
        sum += value
        k += 2

    sum = sqrt(2.0 * sum)
    return sum


def get_norm_factorial(k, m):
    fact = 1
    for i in range(1, 2 * m + 1):
        fact *= k + i
    return fact


def get_pro_norm_cv(m, n, c, cv):
    return get_norm_cv(m, n, c, cv, 1)


def get_obl_norm_cv(m, n, c, cv):
    return get_norm_cv(m, n, c, cv, -1)


def get_pro_norm(m, n, c):
    return get_pro_norm_cv(m, n, c, special.pro_cv(m, n, c))


def get_obl_norm(m, n, c):
    return get_obl_norm_cv(m, n, c, special.obl_cv(m, n, c))

# ---- Shortcut for different types of spheroids

def rad1_cv(m, n, c, cv, type, x):
    if type == 1:
        return special.pro_rad1_cv(m, n, c, cv, x)
    elif type == -1:
        return special.obl_rad1_cv(m, n, c, cv, x)


def rad2_cv(m, n, c, cv, type, x):
    if type == 1:
        return special.pro_rad2_cv(m, n, c, cv, x)
    elif type == -1:
        return special.obl_rad2_cv(m, n, c, cv, x)


def ang1_cv(m, n, c, cv, type, x):
    if type == 1:
        return special.pro_ang1_cv(m, n, c, cv, x)
    elif type == -1:
        return special.obl_ang1_cv(m, n, c, cv, x)


def get_cv(m, n, c, type):
    if type == 1:
        return special.pro_cv(m, n, c)
    elif type == -1:
        return special.obl_cv(m, n, c)


def rad3_cv(m, n, c, cv, type, x):
    return [2.0 * rad1_cv(m, n, c, cv, type, x)[0] + 2.0 * 0j * rad2_cv(m, n, c, cv, type, x)[0],
            2.0 * rad1_cv(m, n, c, cv, type, x)[1] + 2.0 * 0j * rad2_cv(m, n, c, cv, type, x)[1]]

# ---- Generation of A matrices

def get_A11(particle, c1, nmax):
    return get_A(particle, c1, c1, 3, nmax).transpose()


def get_A12(particle, c1, c2, nmax):
    return -get_A(particle, c2, c1, 1, nmax).transpose()


def get_A10(particle, c1, nmax):
    return -get_A(particle, c1, c1, 1, nmax).transpose()


def get_A21(particle, c1, nmax):
    return get_B(particle, c1, c1, 3, nmax).transpose()


def get_A20(particle, c1, nmax):
    return -get_B(particle, c1, c1, 1, nmax).transpose()


def get_A22(particle, c1, c2, nmax):
    return -get_C(particle, c2, c1, nmax).transpose()


def get_a_functions(m, n, c, cv, type, rank, particle):
    if (rank == 1):
        return lambda nu: rad1_cv(m, n, c, cv, type, particle.function(nu))[0]\
                          * ang1_cv(m, n, c, cv, type, nu)[0] / get_norm_cv(m, n, c, cv, type)
    elif(rank == 3):
        return lambda nu: rad3_cv(m, n, c, cv, type, particle.function(nu))[0]\
                          * ang1_cv(m, n, c, cv, type, nu)[0] / get_norm_cv(m, n, c, cv, type)


def get_b_functions(m, n, c, cv, type, rank, particle):
    metric1 = get_integral_metric(particle)
    if (rank == 1):
        return lambda nu: (metric_nu(nu,particle) / metric_psi(nu,particle)
                           * rad1_cv(m, n, c, cv, type, particle.function(nu))[1]
                           * ang1_cv(m, n, c, cv, type, nu)[0] / get_norm_cv(m, n, c, cv, type)
                           - particle.derivative(nu) * metric_psi(nu, particle) / metric_nu(nu, particle)
                             * rad1_cv(m, n, c, cv, type, particle.function(nu))[0]
                             * ang1_cv(m, n, c, cv, type, nu)[1] / get_norm_cv(m, n, c, cv, type)
                              ) / metric1(nu)
    elif(rank == 3):
        return lambda nu: (metric_nu(nu,particle) / metric_psi(nu,particle)
                           * rad3_cv(m, n, c, cv, type, particle.function(nu))[1]
                           * ang1_cv(m, n, c, cv, type, nu)[0] / get_norm_cv(m, n, c, cv, type)
                           - particle.derivative(nu) * metric_psi(nu, particle) / metric_nu(nu, particle)
                             * rad3_cv(m, n, c, cv, type, particle.function(nu))[0]
                             * ang1_cv(m, n, c, cv, type, nu)[1] / get_norm_cv(m, n, c, cv, type)
                              ) / metric1(nu)


def get_c_functions(m, n, c, cv, type, rank, particle):
    eps = particle.eps
    delta = lambda nu: - metric_phi(nu, particle)
    i_zi_t = lambda nu: particle.d / 2.0 * (particle.derivative(nu) * nu + particle.function(nu)) / get_integral_metric(
        particle)(nu)
    return lambda nu: 1.0 / eps * get_b_functions(m, n, c, cv, type, 1, particle)(nu) -\
                      (1.0 / eps - 1.0) * i_zi_t(nu) / delta(nu) * get_a_functions(m, n, c, cv, type, 1, particle)(nu)


def get_integral_metric(particle):
    return lambda nu: sqrt(metric_nu(nu, particle) * metric_nu(nu, particle)\
    + particle.derivative(nu) * particle.derivative(nu)\
      * metric_psi(nu, particle) * metric_psi(nu, particle))


def get_A(particle, c2, c1, rank, nmax):
    A = zeros((nmax, nmax))
    type = particle.type
    m = 1
    for i in range(nmax):
        for k in range(nmax):
            l = i + m
            n = k + m
            cv_l = get_cv(m, l, c2, type)
            cv_n = get_cv(m, n, c2, type)
            func = lambda nu: get_a_functions(m, l, c2, cv_l, type, rank, particle)(nu)\
                              * ang1_cv(m, n, c1, cv_n, particle.type, nu)[0] / get_norm_cv(m, n, c1, cv_n, type)
            A[i][k] = quad(func, -1, 1)[0]
    return A


def get_Z(get_z_functions, particle, c2, c1, rank, nmax):
    Z = zeros((nmax, nmax))
    type = particle.type
    m = 1
    for i in range(nmax):
        for k in range(nmax):
            l = i + m
            n = k + m
            cv_l = get_cv(m, l, c2, type)
            cv_n = get_cv(m, n, c2, type)
            metric = lambda nu: metric_phi(nu, particle) * get_integral_metric(particle)(nu)
            func = lambda nu: get_z_functions(m, l, c2, cv_l, type, rank, particle)(nu)\
                              * ang1_cv(m, n, c1, cv_n, particle.type, nu)[0] / get_norm_cv(m, n, c1, cv_n, type)\
                              * metric(nu)
            Z[i][k] = quad(func, -1, 1)[0]
    return Z


def get_C(particle, c2, c1, nmax):
    return get_Z(get_c_functions, particle, c2, c1, 1, nmax)


def get_B(particle, c2, c1, rank, nmax):
    return get_Z(get_b_functions, particle, c2, c1, rank, nmax)

#-------Solve equation and find solution of scattering by SVM

def get_fullA(particle, c1, c2, nmax):
    return bmat([[get_A11(particle, c1, nmax), get_A12(particle, c1, c2, nmax)],
        [get_A21(particle, c1, nmax), get_A22(particle, c1, c2, nmax)]])


def get_fullB(inputWave, particle, c1, nmax):
    return bmat([[get_A10(particle, c1, nmax)], [get_A20(particle, c1, nmax)]]) * get_Bin(inputWave, nmax)


def get_Bin(inputWave, nmax):
    m = 0
    b = zeros((nmax,1))
    for i in range(nmax):
        b[i] = inputWave.getB(m,i)
    return b

def getSolution(particle, inputWave, c1, c2, nmax):
    x = linalg.solve(get_fullA(particle, c1, c2, nmax), get_fullB(inputWave, particle, c1, nmax))
    b_sca = x[0:nmax + 1]
    b_int = x[nmax + 1:]