import scipy.special as special

from numpy.lib.scimath import sqrt

import f_spheroid
import f_radial

# ---- Calculation of norm for spheroidal angular functions ----

optimize_norm={}
def get_norm_cv(m, n, c, cv, type):
    sum = 0
    k = 0
    key = str(m)+' '+str(n)+' '+str(c)+' '+str(cv)
    if key in optimize_norm:
        return optimize_norm[key]

    d = f_spheroid.sdmn(m, n, c, cv, type)
    isEven = (n - m) % 2
    if(isEven == 1):
        k = 1

    while (k < d.shape[0]):
        if (type == 1):
            value = d[k / 2] * d[k / 2] / (2.0 * k + 2.0 * m + 1.0) * get_norm_factorial(k, m)
        else:
            value = d[k / 2] * d[k / 2] / (2.0 * k + 2.0 * m + 1.0) * get_norm_factorial(k, m)
        sum += value
        k += 2

    sum = sqrt(2.0 * sum)
    optimize_norm[key] = sum
    return sum

def get_norm_factorial(k, m):
    fact = 1l
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

def rad1_cv(m, n, c, type, x):
    eps = 10e-8;
    if type == 1:
        result = f_radial.rad_fun(0, m, n, c, x, eps)
        return result[0], result[1]
    elif type == -1:
        result =  f_radial.rad_fun(1, m, n, c, x, eps)
        return result[0], result[1]


def rad2_cv(m, n, c, type, x):
    eps = 10e-8;
    if type == 1:
        result = f_radial.rad_fun(0, m, n, c, x, eps)
        return result[2], result[3]
    elif type == -1:
        result =  f_radial.rad_fun(1, m, n, c, x, eps)
        return result[2], result[3]


def ang1_cv(m, n, c, cv, type, x):
    if type == 1:
        value = special.pro_ang1_cv(m, n, c, cv, x)
    elif type == -1:
        value = special.obl_ang1_cv(m, n, c, cv, x)
    return value / get_norm_cv(m, n, c, cv, type)


def get_cv(m, n, c, type):
    if type == 1:
        return special.pro_cv(m, n, c)
    elif type == -1:
        return special.obl_cv(m, n, c)

#according to 1.41 Komarov, Slavyanov "Spheroidal funtions"
def rad3_cv(m, n, c, type, x):
    return [rad1_cv(m, n, c, type, x)[0] + 1j * rad2_cv(m, n, c, type, x)[0],
            rad1_cv(m, n, c, type, x)[1] + 1j * rad2_cv(m, n, c, type, x)[1]]

def rad_cv(m, n, c, type, rank, x):
    if rank == 1:
        return rad1_cv(m, n, c, type, x)
    elif rank == 3:
        return rad3_cv(m, n, c, type, x)